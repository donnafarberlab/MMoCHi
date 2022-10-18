from functools import partial
from scipy.signal import find_peaks
import scipy.stats as stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool, get_context
from .logger import logg
from . import utils
try:
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(x, **kwargs):
        return x

try:
    from skfda.preprocessing.registration import landmark_registration_warping, invert_warping
    from skfda.representation.grid import FDataGrid
except ImportError:
    pass
    #raise ImportError('Please install skfda using pip install scikit-fda==0.5')
    
def landmark_register_adts(adata, batch_key = 'donor', data_key = 'protein', key_added = 'landmark_protein',show=False,single_peaks=[],marker_bandwidths={},peak_overrides={},**kwargs):
    ''' Batch correction for ADT expression
    
    Performs negative and positive peak alignment on histograms of ADT expression. Currently expects log or arcsin normalized ADTs (e.g. log2(CP1k)). 
    For now, this method appears to require scikit-fda==0.6. 
    This was developed in parallel with ADTnorm, a similar function expecting arcsin transformed ADTs https://doi.org/10.1101/2022.04.29.489989 
    
    Parameters
    ---
    adata: AnnData object with 
    batch_key: str, column in the .obs of the adata, corresponding to batch information
    data_key: str, key of the dataframe in the .obsm of the adata, corresponding to the ADT expression information
    key_added: str, key of the dataframe in the .obsm of the adata to insert landmark registered ADT expression
    show: bool or int, Whether to show plotted intermediates to better reveal peak detection. Note, this disables parallel processing. Integers (up to 3)
                       correspond to increased verbosity level.
    single_peaks: list of str, columns in adata.obsm[data_key] corresponding to ADTs to only align a single peak of. If multiple peaks are found, the 
                               minimum on the histogram will be used. 
    marker_bandwidths: dict, in the format {'marker_1':0.5}, to override bandwidths used for individual markers
    peak_overrides: dict of dicts, in the format { 'batch_1':{'marker_1':[0.1,0.5]} }, to override peak detection for individual markers in individual batches
    **kwargs: other key word arguments are passed to _detect_landmarks
    '''
    default_bandwidth = 0.2
    assert len(adata.obs_names) == len(set(adata.obs_names)), 'Obs names must be made unique'
    assert adata.obsm[data_key].isna().sum().sum() == 0, f'There cannot be missing values in the adata.obsm["{data_key}"]'
    markers = adata.obsm[data_key].columns
    adata.obsm[key_added] = pd.DataFrame(index=adata.obs_names,columns=markers)
    
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    adata.uns[key_added+"_peaks"] = dict()
    adata.uns[data_key+"_peaks"] = dict()
    for batch_mask, batch in zip(tqdm(batch_masks), batches):
        input_arrays = [np.array(adata.obsm[data_key].loc[batch_mask,marker]) for marker in markers]
        single_peak_list = [marker in single_peaks for marker in markers]
        bandwidths = [marker_bandwidths[marker] if marker in marker_bandwidths.keys() else default_bandwidth for marker in markers]
        overrides = []
        for marker in markers: 
            if batch in peak_overrides.keys() and marker in peak_overrides[batch].keys():
                overrides.append(peak_overrides[batch][marker])
            else:
                overrides.append([])
        barcodes = adata[batch_mask].obs_names
        adata.uns[key_added+"_peaks"][batch] = dict()
        adata.uns[data_key+"_peaks"][batch] = dict()
        if show:
            res = []
            for i, (input_array,single_peak,bandwidth,override) in enumerate(zip(input_arrays,single_peak_list,bandwidths,overrides)):
                res.append(_landmark_registration((input_array,single_peak,bandwidth,override), base_barcodes=barcodes,show=show,info=(str(batch)+" "+markers[i]), **kwargs))
        else:
            fun = partial(_landmark_registration, base_barcodes=barcodes,show=show, **kwargs)
            with get_context("spawn").Pool() as pool:
                res = pool.map(fun, zip(input_arrays,single_peak_list,bandwidths,overrides))
        for marker,r in zip(markers,res):
            adata.obsm[key_added].loc[barcodes,marker] = r[2]
            adata.uns[data_key+"_peaks"][batch][marker] = r[1]
            adata.uns[key_added+"_peaks"][batch][marker] = r[0]
    adata.obsm[key_added] = adata.obsm[key_added].astype(float)
    return adata

def _landmark_registration(array_single_peak_bandwidth_override,base_barcodes,show=False,**kwargs):
    ''' kwargs get sent to detect landmarks TODO add info'''
    array = array_single_peak_bandwidth_override[0]
    single_peak = array_single_peak_bandwidth_override[1]
    bandwidth = array_single_peak_bandwidth_override[2]
    override = array_single_peak_bandwidth_override[3]
    if max(array) == 0:
        return array
    x = np.linspace(0,10,1000)
    barcodes = base_barcodes[array>0]
    ray = array[array>0]
    kde = stats.gaussian_kde(ray,bw_method='scott')
    kde.set_bandwidth(bandwidth)
    y = kde.pdf(x)
    fd = FDataGrid(y, grid_points=x)
    
    if override == []:
        peaks = _detect_landmarks(y,single_peak=single_peak,show=show,**kwargs)
    else:
        if not type(override) == list:
            peaks = list(override)
        else:
            peaks = override
        peaks = [int(round((peak+min(ray)/(max(ray)-min(ray))),2)*100) for peak in peaks]
    
        
    if len(peaks)>1:
        location = [1,3]
        peaks = [x[peaks[0]],x[peaks[-1]]]
    elif len(peaks)==1:
        location,peaks = [1],[x[peaks[0]]]
    else:
        location,peaks = [1],[max(y)] # SHOULD NEVER BE TRIGGERED?
    warp = landmark_registration_warping(fd,[sorted(peaks)],location=location)
    warp = invert_warping(warp)
    ray = warp(ray).reshape(-1)
    # This merge technique might be quite slow. Unsure.
    df = pd.DataFrame(array,index=base_barcodes)
    df[df>0] = np.nan
    df = df.fillna(pd.DataFrame(ray,index=barcodes))
    if show:
        plt.show()
    return (location, peaks, np.array(df[0]))

def _detect_landmarks(y,single_peak=False,landmark_threshold = 0.0025, peak_min_height = .002, show=False,min_dist=30,min_width=20,min_secondary_landmark=200,info=None,primary_prominence=0.003,secondary_prominence=0.0001):
    '''min_secondary_landmark is treated like just a suggestion
    TODO add info'''
    if single_peak:
        landmarks, _ = find_peaks(y,prominence=0.003,width=min_width,height=peak_min_height)
        if len(landmarks) == 0:
            landmarks = [np.argmax(y)]
        landmarks = list(landmarks)
        heights = [y[landmark] for landmark in landmarks]
        for landmark, height in zip(landmarks, heights):
            if height*5 < max(heights):
                landmarks.remove(landmark)
                heights.remove(height)
        primary_landmarks = landmarks = [landmarks[np.argmax(heights)]]
        rejects = secondary_neg_landmarks = secondary_landmarks = []
    else:
        primary_landmarks, _ = find_peaks(y,prominence=0.003,width=min_width,height=peak_min_height)#,rel_height = peak_min_height)
        if int(show) > 2:
            logg.print(["primary",primary_landmarks,_])
        if len(primary_landmarks) == 0:
            primary_landmarks = [np.argmax(y)]
        deriv = np.gradient(y)
        secondary_landmarks, _ = find_peaks(deriv,prominence=0.0008,width=min_width)
        if int(show) > 2:
            logg.print(["secondary",secondary_landmarks,_])
        secondary_landmarks = [landmark for landmark in secondary_landmarks if ((landmark > min(min_secondary_landmark,min(primary_landmarks))) and
                                                                  (abs(deriv[landmark])<landmark_threshold) and
                                                                  (deriv[landmark]<0.001) and
                                                                  (max(abs(deriv[:landmark]))>landmark_threshold/10) and
                                                                  (max(abs(deriv[landmark:]))>landmark_threshold/10))]
        logg.debug((primary_landmarks,secondary_landmarks, max(y)))
        if len(primary_landmarks) + len(secondary_landmarks) < 2:
            min_secondary_landmark=min_secondary_landmark/3
        secondary_neg_landmarks, _ = find_peaks(-deriv,prominence=0.0008,width=min_width)
        if int(show) > 2:
            logg.print(["secondary",secondary_neg_landmarks,_])
        secondary_neg_landmarks = [landmark for landmark in secondary_neg_landmarks if ((landmark > min_secondary_landmark) and # Removes if it's too early, overwrite if a peak alr ID  before then
                                                                  (abs(deriv[landmark])<landmark_threshold) and
                                                                  (deriv[landmark]>-0.001) and
                                                                  (max(abs(deriv[:landmark]))>landmark_threshold/10) and
                                                                  (max(abs(deriv[landmark:]))>landmark_threshold/10))]
        if len(secondary_neg_landmarks) > 1:
            secondary_neg_landmarks = [max(secondary_neg_landmarks)]
        landmarks = list(primary_landmarks) + list(secondary_landmarks) + list(secondary_neg_landmarks)
        rejects = []
        landmarks = np.array(list(set(landmarks)))
        landmarks.sort()
        if int(show) > 2:
            logg.print([landmarks,primary_landmarks,secondary_landmarks,secondary_neg_landmarks])
        while len(landmarks) > 1 and min(np.diff(landmarks)) < min_dist:
            l1_loc = np.argmin(np.diff(landmarks))
            l2_loc = np.argmin(np.diff(landmarks))+1
            l1,l2 = landmarks[l1_loc], landmarks[l2_loc]
            if l1 in primary_landmarks and not l2 in primary_landmarks:
                landmarks = np.delete(landmarks,l2_loc)
                rejects.append(l2)
            elif l2 in primary_landmarks and not l1 in primary_landmarks:
                landmarks = np.delete(landmarks,l1_loc)
                rejects.append(l1)
            elif y[l1] > y[l2]:
                landmarks = np.delete(landmarks,l2_loc)
                rejects.append(l2)
            else:
                landmarks = np.delete(landmarks,l1_loc)
                rejects.append(l1)
        primary_landmarks = [landmark for landmark in primary_landmarks if landmark in landmarks]
        secondary_landmarks = [landmark for landmark in secondary_landmarks if landmark in landmarks and not landmark in primary_landmarks]
        if len(primary_landmarks) > 1:
            landmarks=primary_landmarks
        elif len(primary_landmarks)+len(secondary_landmarks) > 1:
            landmarks=primary_landmarks+secondary_landmarks
    if show:
        fig, axs = plt.subplots(1,2,figsize=(6,3))
        axs[0].plot(rejects,[y[landmark] for landmark in rejects],'yo')
        axs[0].plot(secondary_neg_landmarks,[y[landmark] for landmark in secondary_neg_landmarks],'go')
        axs[0].plot(landmarks,[y[landmark] for landmark in landmarks],'rx')
        axs[0].plot(y)
        axs[0].set_xticks([0,200,400,600,800,1000])
        if not info is None:
            axs[0].set_title(info)
    if int(show) > 1 and not single_peak:
        axs[1].hlines(y=[-landmark_threshold,landmark_threshold],xmin=0, xmax=1000,color='orange')
        axs[1].plot(rejects,[deriv[landmark] for landmark in rejects],'yo')
        axs[1].plot(landmarks,[deriv[landmark] for landmark in landmarks],'rx')
        axs[1].plot(deriv)
        axs[1].set_xticks([0,200,400,600,800,1000])
        axs[1].set_title('Deriv')
    if show:
        plt.show()
    plt.close()
    return landmarks

def stacked_density_plots(adata, marker_list, batch='donor', data_key = ['protein','landmark_protein'],data_key_colors=['b','r'],aspect=3, height=.85, save_fig = None,subsample=1,bw_adjust=0):
    '''Code adapted from https://python.plainenglish.io/ridge-plots-with-pythons-seaborn-4de5725881af'''
    if batch is None:
        batch='None'
        df = pd.DataFrame('All', index=adata.obs_names, columns = [batch])
    else:
        df = pd.DataFrame(adata.obs[batch], index=adata.obs_names, columns = [batch])
    if not pd.api.types.is_list_like(marker_list):
        marker_list = [marker_list]
    if not pd.api.types.is_list_like(data_key):
        data_key = [data_key]
    markname_full_list = []
    for marker in marker_list:
        for dkey in data_key:
            data, markname_full = utils.get_data(adata,marker,dkey,return_source=True)
            df[markname_full] = data
            markname_full_list.append(markname_full)
    df = df.melt(id_vars = batch)
    df = df[df.value>0]
    df.sort_values(['variable',batch],inplace=True)
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
    df = df.sample(frac=subsample)
    g = sns.FacetGrid(df, row=batch, hue=batch,palette='tab10', aspect=aspect, height=height,col='variable',
                     col_order = markname_full_list, row_order = sorted(df[batch].unique()))
    if bw_adjust > 0:
        g.map_dataframe(sns.kdeplot, x='value', fill=True, alpha=.5,gridsize=100, lw=2,bw_adjust=bw_adjust,common_norm=False)
    else:
        g.map_dataframe(sns.histplot, x='value', stat = 'density', bins=40,common_bins=True,alpha=.5)
    g.fig.subplots_adjust(hspace=-.75)
    g.set_titles('')
    g.set(yticks=[], ylabel=None, xlabel=f'Expression')
    g.despine(left=True)
    left_axes = g.axes.flatten('F')[:len(df[batch].unique())]
    for left_ax, batch_id in zip(left_axes, sorted(df[batch].unique())):
        plt.text(x=0,y=0.1,s=batch_id, transform=left_ax.transAxes,ha='right')
    for i, markname_full in enumerate(markname_full_list):
        color = data_key_colors[data_key.index(markname_full.split("_mod_")[-1])]
        plt.text(x=1/(len(markname_full_list))*(i) + .5/(len(markname_full_list)),y=height,s="\n".join(markname_full.replace('__','_').replace('__','_').replace('_mod_','__').split("_")), transform=g.fig.transFigure,ha='center',c=color)
    for i, dkey in enumerate(data_key):
        if (dkey+'_peaks') in adata.uns.keys():
            for j, marker in enumerate(marker_list):
                for k,b in enumerate(sorted(df[batch].unique())):
                    if b in adata.uns[dkey+'_peaks'].keys() and marker in adata.uns[dkey+'_peaks'][b].keys():
                        # g.axes[k][j*len(data_key)+i].text(0,0,(b[0:4],marker,dkey[0]))
                        for peak in adata.uns[dkey+'_peaks'][b][marker]:
                            g.axes[k][j*len(data_key)+i].axvline(peak,ymax=0.1,color='black')
    if not save_fig is None:
        g.savefig(save_fig)
    plt.show()
    return

def density_plot(adata, marker, batch, batch_id, data_key,bw_adjust=0,step=0.1):
    
    data, markname_full = utils.get_data(adata,marker,data_key,return_source=True)
    data = data[adata.obs[batch]==batch_id]
    if bw_adjust > 0:
        ax = sns.kdeplot(x=data, fill=True,gridsize=100, lw=2,bw_adjust=bw_adjust)
    else:
        ax = sns.histplot(x=data, stat = 'density', bins=80)
    for peak in adata.uns[data_key+'_peaks'][batch_id][marker]:
        ax.axvline(peak,ymax=0.1,color='black')
    ax.set_xticks(np.arange(0, round(max(data)+step,1), step),minor=True)
    ax.tick_params(axis='x',bottom=True,which='both')
    ax.set_title(markname_full)
    plt.show()
    return