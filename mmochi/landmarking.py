from tqdm.notebook import tqdm
from tqdm.contrib.concurrent import process_map
from functools import partial

from skfda.preprocessing.registration import landmark_registration_warping, invert_warping
from skfda.representation.grid import FDataGrid
from scipy.signal import find_peaks
import scipy.stats as stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .logger import logg
from . import utils

default_single_peaks = ['Armenian_Hamster_IgG_Isotype_Ctrl', 'Mouse_IgG1_k_isotype_Ctrl', 'Mouse_IgG2a_k_isotype_Ctrl', 'Mouse_IgG2b_k_isotype_Ctrl', 'Rat_IgG1_k_Isotype_Ctrl',
                        'Rat_IgG2a_k_Isotype_Ctrl', 'Rat_IgG2b_k_Isotype_Ctrl', 'anti_human_CD105', 'anti_human_CD119_IFN-g_R_a_chain', 'anti_human_CD11a', 'anti_human_CD112_Nectin-2', 
                        'anti_human_CD122_IL-2Rb', 'anti_human_CD134_OX40', 'anti_human_CD137_4-1BB', 'anti_human_CD138_Syndecan-1','anti_human_CD14','anti_human_CD152_CTLA-4',
                        'anti_human_CD154','anti_human_CD18','anti_human_CD183_CXCR3','anti_human_CD185_CXCR5','anti_human_CD186_CXCR6','anti_human_CD194_CCR4','anti_human_CD196_CCR6',
                        'anti_human_CD1d', 'anti_human_CD209_DC-SIGN','anti_human_CD223_LAG-3', 'anti_human_CD226_DNAM-1','anti_human_CD23','anti_human_CD26','anti_human_CD267_TACI',
                        'anti_human_CD27', 'anti_human_CD270_HVEM_TR2','anti_human_CD272_BTLA', 'anti_human_CD29', 'anti_human_CD303_BDCA-2','anti_human_CD314_NKG2D', 'anti_human_CD319_CRACC', 
                        'anti_human_CD33', 'anti_human_CD335_NKp46', 'anti_human_CD352_NTB-A','anti_human_CD328_Siglec-7','anti_human_CD370_CLEC9A_DNGR1','anti_human_CD38','anti_human_CD39', 
                        'anti_human_CD40','anti_human_CD48','anti_human_CD41', 'anti_human_CD42b','anti_human_CD45','anti_human_CD49b', 'anti_human_CD49d', 'anti_human_CD52', 'anti_human_CD54',
                        'anti_human_CD58_LFA-3','anti_human_CD71','anti_human_CD79b_Igb','anti_human_CD80','anti_human_CD81_TAPA-1','anti_human_CD82','anti_human_CD83', 'anti_human_CD85j_ILT2',
                        'anti_human_CD86','anti_human_CD95_Fas','anti_human_CXCR1','anti_human_HLA-ABC','anti_human_Ig_light_chain_l','anti_human_LOC-1','anti_human_TCR_g_d',
                        'anti_human_TIGIT_VSTIM3', 'anti_mouse_human_CD44', 'anti_human_HLA-DR', 'anti_human_CD107a_LAMP-1','anti_human_CD27','anti_human_CD169_Sialoadhesin_Siglec-1',
                        'anti_human_CD155_PVR','anti_human_CD101_BB27']

def landmark_register_adts(adata, batch_key = 'donor', data_key = 'protein', key_added = 'landmark_protein',show=False,single_peaks=[],marker_bandwidths={},**kwargs):
    ''' Currently expect log2(CP1k) normalized ADTs. I am unsure how important that is
    For now, appears to require scikit-fda==0.6
    TODO add info and cite paper
    
    '''
    default_bandwidth = 0.2
    assert len(adata.obs_names) == len(set(adata.obs_names)), 'Obs names must be made unique'
    assert adata.obsm[data_key].isna().sum().sum() == 0, f'There cannot be missing values in the adata.obsm["{data_key}"]'
    markers = adata.obsm[data_key].columns
    adata.obsm[key_added] = pd.DataFrame(index=adata.obs_names,columns=markers)
    
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    
    for batch_mask, batch in zip(batch_masks, tqdm(batches,total=len(batches))):
        input_arrays = [np.array(adata.obsm[data_key].loc[batch_mask,marker]) for marker in markers]
        single_peak_list = [marker in single_peaks for marker in markers]
        bandwidths = [marker_bandwidths[marker] if marker in marker_bandwidths.keys() else default_bandwidth for marker in markers]
        barcodes = adata[batch_mask].obs_names
        if show:
            arrays = []
            for i, (input_array,single_peak,bandwidth) in enumerate(zip(input_arrays,single_peak_list,bandwidths)):
                arrays.append(_landmark_registration((input_array,single_peak,bandwidth), base_barcodes=barcodes,show=show,info=(str(batch)+" "+markers[i]), **kwargs))
        else:
            fun = partial(_landmark_registration, base_barcodes=barcodes,show=show, **kwargs)
            arrays = process_map(fun, zip(input_arrays,single_peak_list,bandwidths),total=len(input_arrays))
        for marker,array in zip(markers,arrays):
            adata.obsm[key_added].loc[barcodes,marker] = array
    adata.obsm[key_added] = adata.obsm[key_added].astype(float)
    return adata

def _landmark_registration(array_single_peak_bandwidth,base_barcodes,show=False,**kwargs):
    ''' kwargs get sent to detect landmarks TODO add info'''
    array = array_single_peak_bandwidth[0]
    single_peak = array_single_peak_bandwidth[1]
    bandwidth = array_single_peak_bandwidth[2]
    if max(array) == 0:
        return array
    x = np.linspace(0,10,1000)
    barcodes = base_barcodes[array>0]
    ray = array[array>0]
    kde = stats.gaussian_kde(ray,bw_method='scott')
    kde.set_bandwidth(bandwidth)
    y = kde.pdf(x)
    peaks = _detect_landmarks(y,single_peak=single_peak,show=show,**kwargs)
    fd = FDataGrid(y, grid_points=x)
    if len(peaks)>1:
        if y[peaks[0]]*10>y[peaks[-1]]: # If it's a "rare" negative population
            location = [1,3]
        else:
            location = [0.1,1]
        peaks = [x[peaks[0]],x[peaks[-1]]]
    elif len(peaks)==1:
        location,peaks = [1],[x[peaks[0]]]
    else:
        location,peaks = [1],[max(y)]
    warp = landmark_registration_warping(fd,[sorted(peaks)],location=location)
    warp = invert_warping(warp)
    ray = warp(ray).reshape(-1)
    # This merge technique might be quite slow. Unsure.
    df = pd.DataFrame(array,index=base_barcodes)
    df[df>0] = np.nan
    df = df.fillna(pd.DataFrame(ray,index=barcodes))
    return np.array(df[0])

def _test_landmarks(adata,marker, batch, data_key, batch_key,single_peak=False):
    ''' TODO add info'''
    x = np.linspace(0,10,1000)
    ray = adata.obsm[data_key][marker]
    if batch_key is None:
        batches = ['None']
        ray_list = [ray]
    elif batch is None:
        batches=adata.obs[batch_key].unique()
        ray_list = [ray[adata.obs[batch_key]==batch] for batch in batches]
    else:
        batches=[batch]
        ray_list = [ray[adata.obs[batch_key]==batch]]
    for ray,batch in zip(ray_list,batches):
        ray = ray[ray>0]
        kde = stats.gaussian_kde(ray,bw_method='scott')
        kde.set_bandwidth(0.2)
        y = kde.pdf(x)
        peaks = _detect_landmarks(y,show=3,info=batch+' '+marker,single_peak=single_peak)
    return

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
        primary_landmarks = landmarks = [max(landmarks)]
        rejects = secondary_neg_landmarks = secondary_landmarks = []
    else:
        primary_landmarks, _ = find_peaks(y,prominence=0.003,width=min_width,height=peak_min_height)#,rel_height = peak_min_height)
        if int(show) > 2:
            print("primary",primary_landmarks,_)
        if len(primary_landmarks) == 0:
            primary_landmarks = [np.argmax(y)]
        deriv = np.gradient(y)
        secondary_landmarks, _ = find_peaks(deriv,prominence=0.0008,width=min_width)
        if int(show) > 2:
            print("secondary",secondary_landmarks,_)
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
            print("secondary_neg",secondary_neg_landmarks,_)
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
            print(landmarks,primary_landmarks,secondary_landmarks,secondary_neg_landmarks)
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
        fig.show()
    return landmarks