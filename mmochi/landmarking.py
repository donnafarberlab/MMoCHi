from functools import partial
from scipy.signal import find_peaks
import scipy.stats as stats
import numpy as np
import anndata
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from multiprocessing import Pool, get_context
from .logger import logg
from . import utils
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable, List
try:
    from tqdm.auto import tqdm
except ImportError:
    def tqdm(x, **kwargs):
        return x

import os, sys

class _HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        
with _HiddenPrints(): # Attempt to hide an import warning of cython failing, which does not seem to affect performance
    # https://github.com/SheffieldML/GPy/issues/902
    try:
        from skfda.preprocessing.registration import landmark_registration_warping, invert_warping
        from skfda.representation.grid import FDataGrid
    except ImportError:
        try:
            from skfda.preprocessing.registration import landmark_elastic_registration_warping, invert_warping
            from skfda.representation.grid import FDataGrid
        except ImportError:
            pass
            # No longer required to import the package, but is needed to run landmarking functions. Will need to add checking for that
    
def landmark_register_adts(adata: anndata.AnnData, batch_key: str=utils.BATCH_KEY,
                           data_key: str='protein', key_added: str='landmark_protein',
                           show: Union[bool, int]=False, single_peaks:List[str]=[],
                           marker_bandwidths: dict={}, peak_overrides: Union[dict,str]={},
                           inclusion_mask: Optional[Union[List, str]]=None, **kwargs) -> anndata.AnnData:
    ''' Batch correction for expression of all ADTs. 
    
    Performs negative and positive peak alignment on histograms of ADT expression. Currently expects log or arcsine normalized ADTs (e.g. log2(CP1k)). 
    This was developed in parallel with ADTnorm, a similar function which uses arcsine transformed ADTs https://doi.org/10.1101/2022.04.29.489989 
    
    Parameters
    ---
    adata
        object with [batch_key] in .obs, [data_key] log or arcsine normalized data in .obsm to be batch corrected
    batch_key
        column in the .obs of the adata, corresponding to batch information
    data_key
        key of the DataFrame in the .obsm of the adata, corresponding to the ADT expression information
    key_added
        key of the DataFrame in the .obsm of the adata to insert landmark registered ADT expression
    show
        Whether to show plotted intermediates to better reveal peak detection. Note, this disables parallel processing. 
        Integers (up to 3) correspond to increased verbosity level.
    single_peaks
        Columns in adata.obsm[data_key] corresponding to ADTs to only align a single peak of. If multiple peaks are found, the maximum on the histogram will be used. 
    marker_bandwidths
        In the format {'marker_1':0.5}, to override bandwidths used for individual markers
    peak_overrides
        In the format { 'batch_1':{'marker_1':[0.1,0.5]} }, to override peak detection for individual markers in individual batches. 
        To force a positive peak, pass [None,float] as the override. To force a negative peak, you can pass simply [float]. 
        Can also be passed as a str, corresponding to the file path of a JSON saved with this same format.
    inclusion_mask
        Mask events included in landmark registration warping calculation or a column in the .obs with such a mask
    **kwargs 
        Other key word arguments are passed to _detect_landmarks
        
    Returns
    -------
    Adata with batch corrected [data_key], correcting over [batch_key]
    
    '''
    if type(peak_overrides) is str:
        peak_overrides = load_peak_overrides(peak_overrides)
        
    if inclusion_mask is None:
        inclusion_mask = np.array([True] * len(adata))
    if isinstance(inclusion_mask, str):
        inclusion_mask = adata.obs[inclusion_mask]
    
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
                res.append(_landmark_registration((input_array,single_peak,bandwidth,override), base_barcodes=barcodes,show=show,info=(str(batch)+" "+markers[i]),inclusion_mask=inclusion_mask[batch_mask], **kwargs))
        else:
            fun = partial(_landmark_registration, base_barcodes=barcodes,show=show,inclusion_mask=inclusion_mask[batch_mask],  **kwargs)
            with get_context("spawn").Pool() as pool:
                res = pool.map(fun, zip(input_arrays,single_peak_list,bandwidths,overrides))
        for marker,r in zip(markers,res):
            adata.obsm[key_added].loc[barcodes,marker] = r[2]
            adata.uns[data_key+"_peaks"][batch][marker] = r[1]
            adata.uns[key_added+"_peaks"][batch][marker] = r[0]
    adata.obsm[key_added] = adata.obsm[key_added].astype(float)
    return adata

def update_peak_overrides(batch: str, marker: str,
                          update_lower: Optional[Union[float,bool]], 
                          update_upper: Optional[Union[float,bool]],
                          peak_overrides: Union[dict,str]={}, current_peaks: dict=None) -> dict:
    ''' 
    Update peak overrides object for a single batch, marker. Can load it from a JSON if provided a string. 
    To use update_upper and update_lower = False, must provide current_peaks dict (often found in adata.uns)
    
    Parameters
    ---
    batch
        The ID corresponding with the batch of interest
    marker
        The name of the marker to edit overrides on
    update_lower
    update_upper
        If update_lower or update_upper is False, keep current value for lower or upper, respectively. 
        Otherwise, if they are floats, update to that value.
        If update_lower is None, then put None as the lower peak, forcing a single positive. 
        If update_upper is None, then put lower as the only item in the override, forcing a single negative.
        
    Returns
    ---
        peak_overrides: dict
            updated peaks_overrides, which can be passed into update_landmark_register, or saved to JSON using save_peak_overrides    
    '''    
    if type(peak_overrides) is str:
        peak_overrides = load_peak_overrides(peak_overrides)    
    if not type(peak_overrides) is dict:
        peak_overrides = dict()
    if not batch in peak_overrides:
        peak_overrides[batch] = dict()
        
    if update_lower == False:
        update_lower = current_peaks[batch][marker][0]
    if update_upper == False:
        update_upper = current_peaks[batch][marker][-1]
    assert update_lower is None or update_upper is None or (update_upper > update_lower), 'Check values.' 
    assert not update_lower is None or not update_upper is None, 'Cannot both be None.'
    if update_upper is None:
        peak_overrides[batch][marker] = [update_lower]
    else:
        peak_overrides[batch][marker] = [update_lower, update_upper]
        
    return peak_overrides
    
    
                               
def update_landmark_register(adata: anndata.AnnData, batch: str,
                             marker: str, override: Union[List[float],dict,str]={}, 
                             batch_key: str=utils.BATCH_KEY,
                             data_key: str='protein', key_added: str='landmark_protein', 
                             show: Union[bool, int]=False,
                             single_peaks: Union[bool,List[str]]=[],
                             bandwidth: Union[float,dict]=0.2, inclusion_mask: Optional[Union[List, str]]=None,
                             **kwargs)-> anndata.AnnData:
    '''  
    Landmark registration batch correction for ADT expression for a single marker on a single batch. landmark_register_adts() must be run before this. See landmark_register_adts for more details.
    
    Parameters
    ---
    adata
        object with [batch_key] in .obs, [data_key] log or arcsine normalized data in .obsm to be batch corrected
    batch
        Batch to in .obs[batch_key] to define batch to landmark register
    marker
        Corresponding marker of interest
    override
        Can be a listlike of two floats corresponding to the positive and negative peaks, or a dict with an item in { batch:{marker:[0.1,0.5]} }. 
        To force a positive peak, pass [None,float] as the override. To force a negative peak, you can pass simply [float]. 
        Can also be passed as a str, corresponding to the file path of a JSON saved with this dict format.
    batch_key
        column in the .obs of the adata, corresponding to batch information
    data_key
        key of the DataFrame in the .obsm of the adata, corresponding to the ADT expression information
    key_added
        key of the DataFrame in the .obsm of the adata to insert landmark registered ADT expression
    show
        Whether to show plotted intermediates to better reveal peak detection. Note, this disables parallel processing. 
        Integers (up to 3) correspond to increased verbosity level.
    single_peaks
        Columns in adata.obsm[data_key] corresponding to ADTs to only align a single peak of. Can also just be True or False.
    marker_bandwidths
        In the format {marker:0.5}, to override bandwidths used for individual markers, or just a bandwidth number to use that.
    inclusion_mask
        Mask of events included in landmark registration warping calculation or a column in the .obs with such a mask
    **kwargs 
        Other key word arguments are passed to _detect_landmarks
        
    Returns
    -------
    Adata with batch corrected [data_key] over batch_key
    '''
    assert len(adata.obs_names) == len(set(adata.obs_names)), 'Obs names must be made unique'
    assert data_key in adata.obsm.keys(), f'adata.obsm["{data_key}"] not found.'
    assert key_added in adata.obsm.keys() and (key_added+"_peaks" in adata.uns.keys()) and (data_key+"_peaks" in adata.uns.keys()), 'Please run landmark_register_adts first.'
    assert adata.obsm[data_key].isna().sum().sum() == 0, f'There cannot be missing values in the adata.obsm["{data_key}"]'
    assert marker in adata.obsm[data_key].columns, 'Marker not found'
    
    barcodes = adata.obs_names[adata.obs[batch_key] == batch]
    assert len(barcodes) > 1, f'Batch of {batch} in adata.obs["{batch_key}"] not found.'
    
    if inclusion_mask is None:
        inclusion_mask = np.array([True] * len(adata))
    if isinstance(inclusion_mask, str):
        inclusion_mask = adata.obs[inclusion_mask]
    if type(override) is str:
        override = load_peak_overrides(override)   
    if type(override) is dict:
        if batch in override and marker in override[batch]:
            override = override[batch][marker]
        else:
            override = []
    if type(single_peaks) is list:
        single_peaks = marker in single_peaks
    if type(bandwidth) is dict:
        bandwidth = bandwidth[marker]

    inclusion_mask = inclusion_mask[adata.obsm[data_key].index.isin(barcodes)]
    input_array = np.array(adata.obsm[data_key].loc[barcodes,marker])
    res = _landmark_registration((input_array,single_peaks,bandwidth,override), base_barcodes=barcodes,show=show,info=(str(batch)+" "+marker),inclusion_mask=inclusion_mask, **kwargs)
    adata.obsm[key_added].loc[barcodes,marker] = res[2]
    adata.uns[data_key+"_peaks"][batch][marker] = res[1]
    adata.uns[key_added+"_peaks"][batch][marker] = res[0]
    adata.obsm[key_added] = adata.obsm[key_added].astype(float)
    
    return adata


def _landmark_registration(array_single_peak_bandwidth_override: Tuple[np.ndarray, bool, Union[float, int], List[Union[float,int]]],
                           base_barcodes: List[str], show: Union[bool,int]=False, inclusion_mask: Optional[List[bool]]=None,
                           **kwargs) -> Tuple[np.array, np.array, np.array]:
    '''
    Aligns data positive and negative peaks to specified locations
    
    Performs negative and positive peak alignment for one batch of ADT expression. Currently expects log or arcsine normalized ADTs (e.g. log2(CP1k)). 
    This was developed in parallel with ADTnorm, a similar function expecting arcsine transformed ADTs https://doi.org/10.1101/2022.04.29.489989 
    
    Parameters
    ----------
    array_single_peak_bandwidth_override
        Tuple in the order of:
            input_array: numpy array
                Array of expression values for one batch
            single_peak: bool
                If true only look for one peak to correct on
            bandwidths: float
                Mernel density estimation factor to be used in Gaussian KDE
            overrides: list of floats
                Manual peaks location to be used for this batch
    base_barcodes
        DNA barcodes, index aligned to input_array
    show
        Whether to show plotted intermediates to better reveal peak detection. Integers (up to 3) correspond to increased verbosity level.
    inclusion_mask
        Mask of values to use for calculation of peaks (masked adata passed to _detect_landmarks)
    **kwargs
        Other keyword arguments are passed to _detect_landmarks
        
    Returns
    -------
    location:
        Array of position of (negative and) positive population delimiter values
    peaks:
        Array of position of (negative and) positive population delimiters, prior to transformation
    out_arr:
        np.arr with data realigned to have (negative and) positive populations align with [location] values
    '''
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
        if not inclusion_mask is None:
            inclusion_mask = inclusion_mask[array>0]
            kde_masked = stats.gaussian_kde(ray[inclusion_mask],bw_method='scott')
            kde_masked.set_bandwidth(bandwidth)
            y_masked = kde_masked.pdf(x)
        else:
            y_masked = y
        peaks = _detect_landmarks(y_masked,single_peak=single_peak,show=show,**kwargs)
    else:
        if not type(override) == list:
            peaks = list(override)
        else:
            peaks = override
        # peaks must be converted to locations in the array, rather than raw values on the histogram.
        if (peaks[0] is None or peaks[0] == 'None') and not (peaks[-1] is None or peaks[-1] == 'None'):
            peaks = [None,int(round(peaks[-1],2)*100)]
        else:
            peaks = [int(round(peak,2)*100) for peak in peaks] 
        
    if len(peaks)>1:
        if (peaks[0] is None or peaks[0] =='None'): # This allows for an override to be provided for single-positive peaks 
            location = [3]
            peaks = [x[peaks[-1]]]
        else:
            location = [1,3]
            peaks = [x[peaks[0]],x[peaks[-1]]]
    elif len(peaks)==1:
        location,peaks = [1],[x[peaks[0]]]
    else:
        logg.warn('There were no peaks detected. Warning. Setting peak to max of array.')
        location,peaks = [1],[max(y)]
    try:
        warp = landmark_registration_warping(fd,[sorted(peaks)],location=location)
    except NameError:
        warp = landmark_elastic_registration_warping(fd,[sorted(peaks)],location=location)
    warp = invert_warping(warp)
    ray = warp(ray).reshape(-1)
    # This merge technique might be inefficient. Unsure.
    df = pd.DataFrame(array,index=base_barcodes)
    df[df>0] = np.nan
    df = df.fillna(pd.DataFrame(ray,index=barcodes))
    if show:
        plt.show()
    return (location, peaks, np.array(df[0]))

def _detect_landmarks(y: np.array, single_peak: bool=False,
                      landmark_threshold: float=0.0025, peak_min_height: float=.002,
                      show: Union[bool,int]=False, min_dist: Union[int,float]=30,
                      min_width: Union[int,float]=20, min_secondary_landmark: Union[int,float]=200, 
                      info: str=None, primary_prominence: Union[int,float]=0.003, 
                      secondary_prominence: Union[int,float]=0.0008) -> List[Union[int, float]]:
    '''
    Finds the location of peak(s) in the probability density function data 
    
    Uses on scipy.signal.find_peaks to find the locations of landmarks in the inputted data in order to separate (negative and) positive populations.
    
    Parameters
    ----------
    y
        Binned array of values in which to find one or more data peaks 
    single_peak
        If true only try to detect a single peak in the data
    landmark_threshold number:
        Maximum derivative and minimum negative derivative where a peak will be counted as a peak
    peak_min_height
        Minimum height that a peak must be in order to be detected
    show
        Whether to show plotted data to better reveal peak detection. Integers (up to 3) correspond to increased verbosity level.
    min_dist
        Smallest permitted horizontal distance between two peaks
    min_width
        Smallest permitted width of a peak for it to be detected
    min_secondary_landmark
        Smallest permitted value where positive peak can be placed (most used as a suggestion)
    primary_prominence
        Required minimum prominence of peaks. More info at scipy.signal.find_peaks
    secondary_prominence
        A lower minimum prominence to try to find peaks for secondary landmarking. More info at scipy.signal.find_peaks
        
    Returns
    -------
    landmarks
        List of int or float corresponding to (negative and) positive peak location(s)
    '''
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
        secondary_neg_landmarks = [landmark for landmark in secondary_neg_landmarks if ((landmark > min_secondary_landmark) and # Removes if it's too early, 
                                                                  (abs(deriv[landmark])<landmark_threshold) and        # overwrite if a peak alr ID before then
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
        plt.tight_layout()
    plt.close()
    return landmarks

def stacked_density_plots(adata: anndata.AnnData, marker_list: Union[pd.DataFrame, list, tuple],
                          batch_key: str=utils.BATCH_KEY, data_key: List[str]=['protein','landmark_protein'],
                          data_key_colors: Union[List[str], List[Tuple[float]]]=['b','r'],
                          aspect: Union[float,int]=3, height: Union[float,int]=.85, 
                          save_fig: str=None, subsample: Union[float, int]=1,
                          bw_adjust: Union[float,int]=0,exclude_zeroes: bool=True):
    '''
    Method to plot multiple density plots of positive and negative peaks for batches and [data_keys]s with properly placed labels. Method will create density plots for all items in a batch to help visualize the process of realignment (landmark registration)
    
    Code adapted from https://python.plainenglish.io/ridge-plots-with-pythons-seaborn-4de5725881af
    
    Parameters
    ----------
    adata
        AnnData containing [batch] labels in obs, [data_key] in obsm, and contains all items in [markers_list] to plot found using mmc.utils.get_data
    marker_list
        List of markers to be plotted for each [batch], uses mmc.utils.get_data to search for marker info
    batch_key
        Chosen labels to compare. Each item in batch will be its own row in the stacked plot
    data_key
        List of labels to on which to compare [batch]es. Uses mmc.utils.get_data to find labels in data_key in [adata]
    data_key_colors
        Colors to use for text labels for [data_key]s
    aspect
        Aspect value passed to seaborne.FacetGrid function
    height
        Height of the plots. Passed to seaborne.FacetGrid and matplotlib.plt.text functions
    save_fig
        Filepath to save figure to. If none, will not save figure
    subsample
        Fraction of adata data to use for plotting. If less than 1 that fraction will be chosen randomly
    bw_adjust
        Scalar to multiply bandwidth smoothing method used by seaborn.kdeplot. See seaborne.kdeplot for more details
    exclude_zeroes
        If True, only displays non-zero events, which can be useful for visualization of data with many events with zero protein expression
    Returns
    -------
    None. Plots labeled stacked density plots of positive and negative peaks to compare [batch]es across [data_key]s
    
    '''
    if batch_key is None:
        batch_key='None'
        df = pd.DataFrame('All', index=adata.obs_names, columns = [batch_key])
    else:
        df = pd.DataFrame(adata.obs[batch_key], index=adata.obs_names, columns = [batch_key])
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
    df = df.melt(id_vars = batch_key)
    if exclude_zeroes:
        df = df[df.value>0]
    df.sort_values(['variable',batch_key],inplace=True)
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0), 'axes.linewidth':2})
    df = df.sample(frac=subsample)
    g = sns.FacetGrid(df, row=batch_key, hue=batch_key,palette='tab10', aspect=aspect, height=height,col='variable',
                     col_order = markname_full_list, row_order = sorted(df[batch_key].unique()))
    if bw_adjust > 0:
        g.map_dataframe(sns.kdeplot, x='value', fill=True, alpha=.5,gridsize=100, lw=2,bw_adjust=bw_adjust,common_norm=False)
    else:
        g.map_dataframe(sns.histplot, x='value', stat = 'density', bins=40,common_bins=True,alpha=.5)
    g.fig.subplots_adjust(hspace=-.75)
    g.set_titles('')
    g.set(yticks=[], ylabel=None, xlabel=f'Expression')
    g.despine(left=True)
    left_axes = g.axes.flatten('F')[:len(df[batch_key].unique())]
    for left_ax, batch in zip(left_axes, sorted(df[batch_key].unique())):
        plt.text(x=0,y=0.1,s=batch, transform=left_ax.transAxes,ha='right')
    for i, markname_full in enumerate(markname_full_list):
        color = data_key_colors[data_key.index(markname_full.split("_mod_")[-1])]
        plt.text(x=1/(len(markname_full_list))*(i) + .5/(len(markname_full_list)),y=height,
                 s="\n".join(markname_full.replace('__','_').replace('__','_').replace('_mod_','__').split("_")), 
                 transform=g.fig.transFigure,ha='center',c=color)
    for i, dkey in enumerate(data_key):
        if (dkey+'_peaks') in adata.uns.keys():
            for j, marker in enumerate(marker_list):
                for k,b in enumerate(sorted(df[batch_key].unique())):
                    if b in adata.uns[dkey+'_peaks'].keys() and marker in adata.uns[dkey+'_peaks'][b].keys():
                        for peak in adata.uns[dkey+'_peaks'][b][marker]:
                            g.axes[k][j*len(data_key)+i].axvline(peak,ymax=0.1,color='black')
    if not save_fig is None:
        g.savefig(save_fig)
    plt.show()
    return 

def density_plot(adata: anndata.AnnData, marker: str,
                 batch: str, batch_key: str=utils.BATCH_KEY,
                 data_key: str=utils.DATA_KEY,
                 bw_adjust: float=0, step: float=0.1,
                 exclude_zeroes: bool=True):
    '''
    Creates density plot of a single batch for a single marker for a given data_key. Marks location of landmark(s) with small black line(s) on x axis.
    
    Parameters
    ----------
    adata
       AnnData to use with a batch in batch_key in .obs, and data in .obsm[data_key].
    marker
        Marker to be plotted. Uses mmc.utils.get_data to find the marker.
    batch_key
        Category in .obs where batch is.
    batch
        Label in .obs[batch_key] to be plotted.
    data_key
        Label in .obsm to be plotted. 
    bw_adjust
        Scalar to multiply bandwidth smoothing method used by seaborn.kdeplot. See seaborne.kdeplot for more details. 
        Use 0 to plot a histogram.
    step
        Size of x ticks on histogram.   
    exclude_zeroes
        If True, only displays non-zero events, which can be useful for visualization of data with many events with zero protein expression
    '''
  
    data, markname_full = utils.get_data(adata,marker,data_key,return_source=True)
    if exclude_zeroes:
        data = data[(adata.obs[batch_key]==batch) & (adata.obsm[data_key][marker]>0)]
    else:
        data = data[(adata.obs[batch_key]==batch)]
    
    if bw_adjust > 0:
        ax = sns.kdeplot(x=data,fill=True,gridsize=100,lw=2,bw_adjust=bw_adjust)
    else:
        ax = sns.histplot(x=data, stat = 'density', bins=80)
    for peak in adata.uns[data_key+'_peaks'][batch][marker]:
        ax.axvline(peak,ymax=0.1,color='black')
    ax.set_xticks(np.arange(0, round(max(data)+step,1), step),minor=True)
    ax.tick_params(axis='x',bottom=True,which='both')
    ax.set_title(markname_full)
    plt.show()
    return


def density_plot_total(adata: anndata.AnnData, marker: str,
                       batch: str, batch_key: str=utils.BATCH_KEY,
                       data_key: str=utils.DATA_KEY, 
                       bw_adjust: float=0, step: float=0.1,
                       weights: Optional[float]=None):
    '''
    Plots density of a single marker on a single batch in front of the density of this marker for the whole dataset.
    Creates density plot of a single batch for a single marker for a given data_key. Marks location of landmark(s) with small black line(s) on x axis.
    Plots this in front of the density plot for the entire dataset. This is helpful when adding a single batch to a larger, pre-landmarked dataset. 
    
    Parameters
    ----------
    adata
       AnnData to use with a batch in batch_key in .obs, and data in .obsm[data_key].
    marker
        Marker to be plotted. Uses mmc.utils.get_data to find the marker.
    batch_key
        Category in .obs where batch is.
    batch
        Label in .obs[batch_key] to be plotted.
    data_key
        Label in .obsm to be plotted. 
    bw_adjust
        Scalar to multiply bandwidth smoothing method used by seaborn.kdeplot. See seaborne.kdeplot for more details.
    step
        Size of x ticks on histogram.    
    weight
        See seaborn.histplot documentation, the weights for the overall dataset to manually normalize the histogram
        heights of the batch and the total dataset. Passing None calculates n_batch_events/n_total_events
    '''
    data, markname_full = utils.get_data(adata,marker,data_key,return_source=True)
    data_max = max(data)
    data_batch = data[adata.obs[batch_key]==batch]
    data = data[data>0]
    data_batch = data_batch[data_batch>0]
    if weights is None:
        weights = len(data_batch)/len(data)
    fig, ax = plt.subplots()
    if bw_adjust > 0:
        sns.histplot(x=data,  bins=80, alpha=0.5, color='grey', weights=weights)
        sns.histplot(x=data_batch,  bins=80, fill=True, alpha=0.5, kde=True, color="#CB181D", kde_kws={'bw_adjust':bw_adjust})
    else:
        sns.histplot(x=data,  bins=80, alpha=0.5, color='grey', weights=weights)
        sns.histplot(x=data_batch,  bins=80, alpha=0.5, kde=False, color="#CB181D")
    for peak in adata.uns[data_key+'_peaks'][batch][marker]:
        ax.axvline(peak,ymax=0.1,color='black')
    ax.set_xticks(np.arange(0, round(data_max+step,1), step),minor=True)
    ax.tick_params(axis='x',bottom=True, which='both')
    ax.set_title(markname_full)
    plt.xlim(0,data_max)
    plt.show()
    return

def save_peak_overrides(path: str, peak_overrides: dict):
    '''
    Saves peak overrides to a JSON file for easy loading.
    
    Parameters
    ----------
    path
        Location to save to, including the file name and .json
    peak_overrides
        Dictionary of peak overrides to save
    '''
    import json
    with open(path, 'w') as fp:
        json.dump(peak_overrides, fp, skipkeys=True, sort_keys=True, indent=4)
    return

def load_peak_overrides(path: str):
    '''
    Loads peak overrides from a JSON file.
    
    Parameters
    ----------
    path
        Location to load file from, including the file name and .json
    
    Returns
    -------
    Dictionary of peak overrides from the JSON file.
    '''    
    import json
    with open(path, "r") as fp:
        peak_overrides = json.load(fp)
    return peak_overrides