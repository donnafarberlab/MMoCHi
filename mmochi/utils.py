import anndata
import sys
import os
import pandas as pd
import numpy as np
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import scanpy as sc
import gc
from .logger import logg
import warnings
import matplotlib.pyplot as plt

DATA_KEY = 'landmark_protein'
BATCH_KEY = 'batch'
MODALITY_COLUMN = None

def convert_10x(adata: anndata.AnnData, drop_totalseq: bool=True, data_key: str='protein') -> anndata.AnnData:
    '''
    Convert default 10X data to an AnnData with protein in the .obsm, and gene expression in the .X
    
    Parameters
    ----------
    adata
        10X provided AnnData object
    drop_totalseq
        If true, removes '_TotalSeq' from the end of .obsm[data_key] values
    data_key
        Name of data in .obsm to convert 
        
    Returns 
    -------
        AnnData object with 'Antibody Capture' data moved to .obsm[data_key]
    '''
    adata.obsm[data_key] = pd.DataFrame(index=adata.obs_names)
    for i in adata.var_names[adata.var.feature_types=='Antibody Capture']:
        adata.obsm[data_key][i] = adata.obs_vector(i)
    adata = adata[:,adata.var.feature_types!='Antibody Capture'].copy()
    if drop_totalseq:
        adata.obsm[data_key].columns = pd.Series(adata.obsm[data_key].columns).str.split('_TotalSeq',expand=True)[0].values
    return adata

def _default_kwargs(kwargs: dict, defaults: dict) -> dict:
    '''
    Helper function to combine lists of kwargs with defaults, if no kwargs is specified will use defaults instead
    '''
    for key in defaults.keys():
        if key not in kwargs.keys():
            kwargs[key] = defaults[key]
    return kwargs

def intersect_features(adatas: List[anndata.AnnData], data_key: str=DATA_KEY) -> List[anndata.AnnData]:
    '''
    Subsets each AnnData object to only genes and values of data_key found in every AnnData objects. Note: multiple valid obsm_keys is not supported.
    
    Parameters
    ----------
    adatas
        Each AnnData object must have gene names in .var_names and must have the specified data_key in .obsm
    data_key
        Each AnnData object's .obsm must have data_key as one of its fields
        
    Returns 
    -------
    adatas
        List of AnnData objects with same fields as inputted object, but containing information on only genes and 'data_key' that exists in every object
    '''
    intersection_GEX = list(set.intersection(*(set(val) for val in [adata.var_names for adata in adatas])))
    if not data_key is None:
        obsm_key = list(set(data_key).intersection(set(adatas[0].obsm.keys())))
        if len(obsm_key) > 1:
            logg.warn(f'Multiple valid data_keys {obsm_key} not supported, will use {obsm_key[0]}')
            obsm_key = obsm_key[0]
        if len(obsm_key) != 0:
            intersection_ADT = list(set.intersection(*(set(val) for val in [adata.obsm[data_key].columns for adata in adatas])))
    for i,adata in enumerate(adatas):
        adata = adata[:,intersection_GEX].copy()
        if not data_key is None and len(obsm_key) != 0:
            adata.obsm[data_key] = adata.obsm[data_key][intersection_ADT].copy()
        adatas[i] = adata
    return adatas

def preprocess_adatas(adatas: Union[anndata.AnnData, List[anndata.AnnData], str, List[str]]=None, 
                      convert_from_10X: bool=True,
                      make_unique: bool=True,
                      drop_totalseq: bool=True, 
                      intersect: bool=True, 
                      backup_urls: Union[str, List[str]]=None,
                      data_key: str='protein',
                      log_CP_GEX: int=1e4,
                      log_CP_ADT: int=1e3) -> List[anndata.AnnData]:
    
    '''
    Function to load and preprocess adatas from either filename(s) or backup_url(s). 
    
    Parameters 
    ----------
    adatas
        Individual or list of filepaths or AnnData objects 
    convert_from_10X
        Determines whether to convert from 10x data style
    make_unique
        Determines whether to make labels in AnnData objects unique by adding '-1' etc. to duplicates
    drop_totalseq
        Determines whether to drop '_TotalSeq' from the data_key labels
    intersect
        Determines whether to subset each DataFrame to 'data_key' and gene expression data only found on all DataFrames
    backup_urls
        url or urls to the AnnData objects to be used if the filepath does not work
    data_key
        Must be in .obsm, name used for data intersection and adt normalization
    log_CP_GEX
        Normalizes Gene expression counts to log(counts/total_counts*log_CP_GEX + 1), default 10000
    log_CP_ADT
        Normalizes ADT counts to log(counts/total_counts*log_CP_ADT + 1), default 1000
        
    Returns
    -------
    List of AnnData object(s)
        if make unique:
            containing unique obs and var labels, with duplicates named '-1', '-2' etc. 
        if convert_from_10x:
            with protein data moved to the .obsm
        if intersect:
            containing only the subset of genes and 'data_key' that exists in all objects
        if log_CP_GEX:
            log normalized gene counts per x value 
        if log_CP_ADT: 
            log normalized 'data_key' counts per x value 
     '''
    warnings.filterwarnings("ignore", category=UserWarning) # Ignores var_names not unique warnings
    if isinstance(adatas, anndata.AnnData) or isinstance(adatas,str):
        adatas = [adatas]
    if isinstance(backup_urls, str):
        backup_urls = [backup_urls]
    if backup_urls is None:
        backup_urls = [None] * len(adatas)
    for i,(adata,backup_url) in enumerate(zip(adatas,backup_urls)):
        if isinstance(adata,str):
            if backup_url is None and adata[-4:] == 'h5ad':
                adata = anndata.read_h5ad(adata)
            else:
                adata = sc.read_10x_h5(adata, backup_url=backup_url, gex_only=False)
                
        if make_unique:
            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            
        if convert_from_10X:
            adata = convert_10x(adata,drop_totalseq)
        
        if drop_totalseq and make_unique and convert_from_10X:
            assert len(adata.obsm[data_key].columns) == len(adata.obsm[data_key].columns.unique()), 'Removing TotalSeq tags made antibodies not unique, try' + \
                                                                                                     'again without `drop_totalseq = True`'
        adatas[i] = adata
        
    if intersect:
        adatas = intersect_features(adatas,data_key)
        adatas = adatas
    for i, adata in enumerate(adatas):
        if log_CP_GEX > 0 and ((adata.X.astype(int) != adata.X).nnz==0):
            adata.layers['counts'] = adata.X.astype(int).copy()
            sc.pp.normalize_total(adata, target_sum = log_CP_GEX)
            sc.pp.log1p(adata)
        if log_CP_ADT > 0 and ((data_key in adata.obsm.keys()) and all(adata.obsm[data_key].astype(int) == adata.obsm[data_key])):
            adata.obsm[data_key+"_counts"] = adata.obsm[data_key].copy()
            adata_prot = anndata.AnnData(X= adata.obsm[data_key])
            sc.pp.normalize_total(adata_prot, target_sum = log_CP_ADT)
            sc.pp.log1p(adata_prot)
            adata.obsm[data_key] = pd.DataFrame(adata_prot.X,columns=adata.obsm[data_key].columns,index=adata.obs_names)
        adatas[i] = adata
    return adatas


def batch_iterator(adata: anndata.AnnData, batch_key: str,
                   sort: bool=True) -> Tuple[List[List[bool]], List[str]]:
    '''Generates a series of masks, for each different batch in batch_key and its corresponding batch name
    
    Parameters
    ----------
    adata
    batch_key
        Column within the adata.obs that delineates batches
    sort
        Whether to sort the batch names according to the sorted() function
    
    Returns 
    -------
    batch_masks
        Boolean masks for each of the batches in batch_key
    batches
        The names of the batches in the batch_key column of adata.obs
    '''
    batch_masks, batches = [], []
    if batch_key is None:
        batch_masks.append(pd.Series([True] * adata.n_obs,index=adata.obs_names))
        batches.append(None)
    else:
        logg.info(f'Running with batch {batch_key}')
        items = sorted(adata.obs[batch_key].unique()) if sort else adata.obs[batch_key].unique()
        for batch_id in items:
            batch_masks.append(adata.obs[batch_key] == batch_id)
            batches.append(batch_id)
    return batch_masks, batches

def obsm_to_X(adata: anndata.AnnData, data_key: str=DATA_KEY) -> anndata.AnnData:
    '''
    Makes the .obsm[data_key] of an AnnData object into its .X, if multiple data_key values provided, will merge first found key.
    
    Parameters
    ----------
    adata
    data_key
        key in the adata.obsm to convert to the .X of this new adata
    Returns 
    -------
    adata_new
        Object with .X of AnnData.obsm[data_key]
    '''
    if isinstance(data_key, list):
        data_key = data_key[0]

    adata_new = anndata.AnnData(adata.obsm[data_key],adata.obs,pd.DataFrame(index=adata.obsm[data_key].columns)) # dtype=float
    return adata_new

def _marker(keyword: str, cols: List[str], allow_multiple: bool=False) -> Union[List[str], str]:
    '''
    Function to search for a marker in a list or listlike of markers
    
    Parameters
    ----------
    keyword
        Name of the marker you wish to search for
    cols
        List of marker names in which to search for [keyword]
    allow_multiples
        Whether to return a list of matches if there are multiple matches or throw an error
        
    Returns 
    -------
    if allow_multiples 
        marker_names
            Names of the marker in cols
    else 
        marker_name
            Name of the marker in cols
    '''
    keyword = keyword.lower()
    marker_name = [i for i in cols if keyword in i.lower()]
    
    assert len(marker_name) > 0, "\"" + keyword + "\" not found in marker list."
    if allow_multiple: # If we're allowing multiple just shove out whatever
        return marker_name
    elif len(marker_name) > 1: # if there're more than one result
        # Look for an exact-exact match
        for i in marker_name: 
            if i.lower() == keyword.lower(): 
                return i
        for i in marker_name:
            if i.lower().split(keyword)[1].split('_')[0] == '': 
                return i
        # If no single exact match is found
        raise ValueError(f" '{keyword}' found multiple times in marker list: {marker_name}")
    else:
        return marker_name[0]

def marker(adata: anndata.AnnData, parameter: Union[str, List[str]], 
           data_key: Optional[Union[str,list]]=DATA_KEY, allow_multiple: bool=False) -> Union[str, List[str]]:
    '''
    Lookup a marker name within the .X, the .obsm[data_key], or .var[MODALITY_COLUMN] == data_key.
    
    Parameters
    ----------
    adata
        Object with data in adata.X, adata.obsm[data_key], or where adata.var[MODALITY_COLUMN] == data_key
    parameter
        Parameters to search for in adata
    data_key
        Location in the .obsm or modality in .var to look for marker name
        If None, searches adata.var_names
    allow_multiple
        Whether to return a list of matches if there are multiple matches or throw an error
        
    Returns 
    -------
        marker_name
            Name(s) of the marker in the adata object
    '''
    
    cols_list = []
    if isinstance(data_key, list) or not data_key is None:
        if not isinstance(data_key, list):
            data_key = [data_key]
        for d_key in data_key:
            try: #.var
                if MODALITY_COLUMN is None or MODALITY_COLUMN not in adata.var.columns:
                    assert False, 'no valid modailty column'
                cols = adata.var_names[adata.var[MODALITY_COLUMN] == d_key]
                assert len(cols) > 0, f'data_key {d_key} not found in .var[{MODALITY_COLUMN}]'
                cols_list.append(cols)
            except AssertionError:
                try:
                    cols_list.append(adata.obsm[d_key].columns.to_list())
                except AssertionError:
                    pass
    else:
        cols_list.append(adata.var_names.to_list())

    if not isinstance(parameter, list): # Run recursively if it's a list of inputs
        parameter = [parameter]
    
    results = []
    for i, word in enumerate(parameter):
        for cols in cols_list:
            try:
                results.append(_marker(word, cols, allow_multiple))
                break
            except AssertionError:
                pass
            
    if allow_multiple:
        return sum(results,[])
    elif len(results)==1:
        return results[0]
    elif len(results)==0:
        assert False, "No markers found"
    else:
        return results
        
def umap_thresh(adata: anndata.AnnData, h,
                markers: Union[str,List[str]]=None, batch_key: str=BATCH_KEY,
                data_key: Optional[Union[list,str]]=DATA_KEY, umap_basis: str='X_umap',
                cmap: Optional[List[str]]=None, external_holdout: bool=False, key_added: str = 'lin',**kwargs):
    '''
    Plots UMAPs for the listed markers with thresholded expression data for the markers overlayed on top. 
    This visualization can be useful to see where protein and genes are expressed in the UMAP.
    
    Parameters
    ----------
    adata
        Object with .obsm['X_umap'] with UMAP representation and thresholds to color the umap
    h
        hierarchy object to draw classification markers from
    markers
        Marker(s) to create UMAPs for
        If None (default), creates UMAP for all markers in h
    batch_key
        Name of batch in adata to use for UMAP, uses batch_iterator to find key
    data_key
        Location in the .obsm and .var[MODALITY_COLUMN] to look for marker name
        If None, searches adata.var_names
    umap_basis
        Passed to scanpy.pl.embedding as basis
    cmap
        List of colors or color map information used to color the feature expression on the UMAP
    external_holdout
        If external hold out was defined in adata.obsm[key_added], removes external hold out from automatic threshold calculations and from graphs used for manual thresholding
    key_added
        If external_holdout is true, the place in adata.obsm[key_added] to search for the external hold out column (bool T F column used to indicate whether an event should be set aside for hold out)
    kwargs
        Arguments to be passed to scanpy.pl.embedding
    
    '''
    if external_holdout:
        try:
             external_holdout_mask = adata.obsm[key_added]['external_holdout']
        except:
            assert False, f'Did not create external hold out in adata.obsm[{key_added}]'
    else:
        external_holdout_mask = [True] * len(adata)

    if markers is None:
        all_markers = []
        for classification in h.get_classifications():
            for mm in h.tree[classification].data.markers:
                if mm in adata.obs.columns:
                    logg.warn(f"{mm} found in adata.obs.columns, skipping...")
                else:
                    all_markers.append(mm)
        markers = list(set(all_markers))
    elif isinstance(markers, str):
        markers = [markers]
    
    for mm in markers:
        if mm.endswith("_lo") or mm.endswith("_hi"):
            m = mm[:-3]
        else:
            m = mm
        to_del = False
        try:
            m = marker(adata,m,data_key)+"_mod_"+data_key
            try: 
                assert MODALITY_COLUMN in adata.var.keys()
                assert sum(adata.var[MODALITY_COLUMN] == [m.split('_mod_')[0]]) > 0
                adata.obs[m] = adata[:,adata.var[MODALITY_COLUMN] == [m.split('_mod_')[0]]]
            except:
                adata.obs[m] = adata.obsm[data_key][m.split('_mod_')[0]]
            to_del = True
        except:
            m = marker(adata,m.split('_gex')[0],None)
        data = get_data(adata[external_holdout_mask],m,data_key)
        for mask,b in zip(*batch_iterator(adata[external_holdout_mask],batch_key)):
            t = h.get_threshold_info(mm,None,b,flexible_batch=True)[1]
            t = sorted(t)
            t[-1] = min(t[-1],max(data))
            t[0] = min(t[0],max(data))
            t[0] = max(t[0],min(data))
            t[-1] = max(t[-1],min(data))
            if not type(cmap) is str:
                from matplotlib.colors import ListedColormap
                n_o = np.ceil(2560*(max(data)-max(t))/max(data))
                n_g = np.ceil(2560*(max(t)-min(t))/max(data))
                n_b = np.ceil(2560*min(t)/max(data))
                try: 
                    from matplotlib import cm
                    orange = cm.get_cmap('Oranges')
                    blue = cm.get_cmap('Blues')
                    green = cm.get_cmap('Greens')
                except:
                    orange = plt.get_cmap('Oranges')
                    blue = plt.get_cmap('Blues')
                    green = plt.get_cmap('Greens')                    
                newcolors =  list(blue(np.linspace(0, 1, int(n_b)))[int(np.ceil(0.3*n_b)):int(np.ceil(0.860*n_b))]) + \
                             list(green(np.linspace(0, 1, int(n_g)))[int(np.ceil(0.3*n_g)):int(np.ceil(0.860*n_g))]) + \
                             list(orange(np.linspace(0, 1, int(n_o)))[int(np.ceil(0.3*n_o)):int(np.ceil(0.860*n_o))])
                newcolors = np.insert(newcolors,0,np.array([.95,.95,.95,1]),axis=0)
                cmap = ListedColormap(newcolors)

            sc.pl.embedding(adata[external_holdout_mask][mask],basis = umap_basis,color=m,cmap=cmap,vmax=max(data),vmin=0,title=f'{b} {mm}',**kwargs)
            plt.show()
            gc.collect()
        if to_del:
            del adata.obs[m]
    return

def umap_interrogate_level(adata: anndata.AnnData, level: str, batch_key: str=None,
                            key_added: str='lin', umap_basis: str='X_umap',
                            cmap: Optional[List[str]] = None, **kwargs):
    '''
    Plots UMAPs showing events selected by high-confidence thresholding and used for training and breaks down annotation confusion onto the UMAP. X_umap must be generated before running this function.
        
    Parameters
    ----------
    adata
        AnnData object to use for plotting
    level
        Level of the classifier to interrogate
    batch_key
        Name of batch in adata to use for UMAP, uses batch_iterator to find key
    key_added
        Key in the adata.obsm that contains the classification results
    umap_basis
        Passed to scanpy.pl.embedding as basis
    cmap
        Color map information used to color 
    kwargs
        Sent to scanpy.pl.embedding
    '''
    batch_masks, batches = batch_iterator(adata,batch_key)
    adata.obs['High-Confidence Thresh'] = adata.obsm[key_added][f'{level}_hc'].astype(str).replace({'nan':None})
    adata.obs['High-Confidence Thresh'] = pd.Categorical(adata.obs['High-Confidence Thresh'], 
                                                         categories = sorted(set(adata.obs['High-Confidence Thresh'])-set([None,'?',np.nan]))+['?'])
    adata.obs['Training Data'] = adata.obs['High-Confidence Thresh'].copy()
    adata.obs.loc[~adata.obsm[key_added][f'{level}_train'],'Training Data'] = None
    adata.obs['Training Counts'] = adata.obsm[key_added][f'{level}_traincounts'].replace({0:np.nan})
    
    adata.obs['HC -> Class'] = (adata.obsm[key_added][f'{level}_hc'].astype(str) +" -> "+ adata.obsm[key_added][level+'_class'].astype(str))
    adata.obs['HC -> Class'] = adata.obs['HC -> Class'].replace({'nan -> nan':None})
    adata.obs['HC -> Class'] = pd.Categorical(adata.obs['HC -> Class'], categories = sorted(set(adata.obs['HC -> Class'])-set([None])))

    for batch_mask, batch in zip(batch_masks, batches):
        sc.pl.embedding(adata[batch_mask],basis = umap_basis,color=['High-Confidence Thresh','Training Data','Training Counts','HC -> Class'],
                        cmap=cmap, ncols=2, **kwargs)
        gc.collect()
    adata.obs = adata.obs.drop(['High-Confidence Thresh','Training Data','Training Counts','HC -> Class'],axis=1)
    return
        
def get_data(adata: anndata.AnnData,
             parameter: str,
             preferred_data_key: Optional[Union[List,str]]=None,
             return_source: bool=False,
             force_preferred: bool=False) -> np.array:
    '''
    Searches an AnnData object along its .var, .var[MODALITY_COLUMN], .obs, .layers, and .obsm[preferred_data_key] for a specified parameter. 
        
    Parameters
    ----------
    adata
        Object to search through     
    parameter
        string to search the AnnData object for. Append:
        "_gex" to limit checking to the .var,
        "_obs" to limit checking to the .obs or where .var[MODALITY_COLUMN] == 'obs',
        "_mod_" followed by an var name in .var[MODALITY_COLUMN] or an obs key.
        
        Note, these symbols overwrite preferred data key, and it is assumed
        that these symbols are not used elsewhere for variable names
    preferred_data_key
        Used to point to key(s) in .layers, .var[MODALITY_COLUMN], or in .obsm to check in
    return_source
        Whether to return where the parameter was found in addition to the list of matches
    force_preferred
        Whether to raise an error if the parameter is not found in the preferred_data_key
    
    Returns
    -------
    data
       Array of data in the AnnData object whose label matches parameter
    '''
            
    data = None
    if "_mod_" in parameter:
        force_preferred = True
        preferred_data_key = parameter.split('_mod_')[1]
        parameter = parameter.split('_mod_')[0]
    if "_gex" in parameter or "_obs" in parameter:
        preferred_data_key = None

    if isinstance(preferred_data_key, list):
        for d_key in preferred_data_key:
            try:
                res = get_data(adata,parameter,d_key,return_source,force_preferred)
            except:
                res = None
            if not res is None:
                return res  
        if res is None:
            if force_preferred:
                raise ValueError(f'Did not find {parameter} for any data_key: {preferred_data_key}') 
                    
                    
    if preferred_data_key is not None:   
        assert int(preferred_data_key in adata.obsm.keys()) + int(preferred_data_key in adata.layers.keys()) + int(MODALITY_COLUMN in adata.var.columns and preferred_data_key in adata.var[MODALITY_COLUMN]) <= 1, \
        f'Found {preferred_data_key} in multiple locations in the adata.obsm, adata.X, and adata.X.layers(), please inspect AnnData object and rectify.'
        
        if MODALITY_COLUMN in adata.var.columns:
            try:
                data = np.array(adata[:,adata.var[MODALITY_COLUMN] == preferred_data_key].obs_vector(parameter))
            except:
                try:
                    parameter = _marker(parameter,set(adata.obsm[MODALITY_COLUMN].values),False)
                    data = np.array(adata[:,adata.var[MODALITY_COLUMN] == preferred_data_key].obs_vector(parameter))
                except:
                    pass   
                
        if data is None:      
            if preferred_data_key in adata.layers.keys(): # Try checking layers and obsm
                try:
                    data = np.array(adata.obs_vector(parameter,layer=preferred_data_key))
                except:
                    parameter = _marker(parameter,adata.var_names, False)
                    data = np.array(adata.obs_vector(parameter,layer=preferred_data_key)) 

            elif preferred_data_key in adata.obsm.keys():
                try:
                    data = np.array(adata.obsm[preferred_data_key][parameter])
                except:
                    try:
                        parameter = _marker(parameter,adata.obsm[preferred_data_key].columns,False)
                        data = np.array(adata.obsm[preferred_data_key][parameter])
                    except:
                        pass

                
    if not data is None:
        if return_source:
            return data, parameter + '_mod_' + preferred_data_key
        else:
            return data
    
    if force_preferred:
        raise ValueError(f'Did not find {preferred_data_key} in adata.var[{MODALITY_COLUMN}], adata.obsm, or adata.layers\n'+
                         f'Or did not find {parameter} in that data')
    else:
        assert parameter not in adata.var_names or parameter not in adata.obs.keys(), \
        'Found {parameter} in both adata.obs and adata.var, please inspect AnnData Object and rectify.'

        try: # Try checking .X
            try:
                if parameter.split("_gex")[0] in adata.var_names:
                    data = np.array(adata.obs_vector(parameter.split("_gex")[0]))
                else:
                    assert False
            except:
                parameter = _marker(parameter.split("_gex")[0],adata.var_names, False)
                data = np.array(adata.obs_vector(parameter))
        except:
            pass
        
        if not data is None:
            if return_source:
                if MODALITY_COLUMN in adata.var.columns:
                    return data, parameter.split("_gex")[0] + '_mod_' + adata.var.loc[parameter.split("_gex")[0],MODALITY_COLUMN]
                else:
                    return data, parameter.split("_gex")[0]+'_gex'
            else:
                return data
        
        try: # Try checking .obs
            try:
                data = np.array(adata.obs[parameter.split("_obs")[0]])
            except:
                parameter = _marker(parameter.split("_obs")[0],adata.obs.columns.to_list(), False)
                data = np.array(adata.obs[parameter])
        except:
            pass
        
        if not data is None:
            if return_source:
                return data, parameter.split("_obs")[0]+'_obs'
            else:
                return data

    if preferred_data_key is None:
        raise ValueError(f'Did not find {parameter} in adata.var or adata.obs') 
    else:
        raise ValueError(f'Could not find {parameter} in adata.')

def _list_tools(list1: list, operator: str,
               list2: list=[]) -> list:
    '''
    Helper functions to combine lists. Takes in two lists, with an operator for how to join the lists
    
    Parameters
    ----------
    list1
        A list
    Operator
        '==' followed by a value - Returns a boolean list of whether each value is equal to the value after == \n
        '+'  - Returns the sum of the lists \n
        '-'  - Returns the values in list1 that are not in list2 \n
        '*'  - Returns the all permutations of sum of one element from list1 and one element from list2, in the order of [list1[0] + list2[0], list1[1] + list2[0], list1[2] +list2[0], ... , list1[n] + list2[m]]
    list2
        Another list
    
    '''
    if '==' in operator:
        return [i == operator.split('==')[1] for i in list1]
    elif operator == '+':
        return list1 + list2
    elif operator == '-':
        return [i for i in list1 if not i in list2]
    elif operator == "*":
        new_list2 = list2*len(list1)
        new_list1 = [item for item in list1 for i in range(len(list2))]
        return [i+j for i,j in zip(new_list1,new_list2)]
    else:
        return None
    
def generate_exclusive_features(adata_list: Union[List[str], List[anndata.AnnData]], 
                                data_key: Optional[Union[str,list]]=DATA_KEY):
    '''
    Reads in adata objects from a list of paths, or takes in a list of adata objects
    and finds the features in all objects
    
    Parameters
    ----------
    adata_list
        List of objects in which to search for common features
    data_key
        Optional key in the adata.obsm to use for features, will be intersected in addition to .var_names
        
    Returns
    -------
    features
        features that appear in all AnnData objects
    '''
    features = []
    for adata in adata_list:
        if isinstance(adata, str):
            adata = anndata.read_h5ad(adata)
        if data_key is None:
            features_one = list(adata.var_names)
        else:
            if not isinstance(data_key, list):
                data_key = [data_key]
            features_one = list(adata.var_names)
            for d_key in data_key:
                try:
                    features_one.extend(list(adata.obsm[data_key].columns))
                except:
                    pass
        if features != []:
            features = [f for f in features_one if f in features]
            logg.print("intersection taken")
        else:
            logg.print("No itersection attempted")
            features = features_one
    return features

def _validate_names(names: List[str], reserved_strings: List[str]=['_obs','_gex','_mod_']):
    '''
    Checks names to ensure none of the reserved strings are contained in them

    Parameters
    ----------
    names
        iterable of strings to check
    reserved_string
        list of substrings that can't be contained in names
    '''
    for name in names:
        for reserved_string in reserved_strings:
            assert reserved_string not in name, f'invalid name {name}, contains {reserved_string}'
    return

