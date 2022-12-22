import anndata
import sys
import os
import pandas as pd
import numpy as np
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import scanpy as sc
from .logger import logg

def convert_10x(adata: anndata.AnnData, drop_totalseq: bool = True, data_key: str = 'protein'):
    '''
    Convert default 10X data to an anndata with protein in the .obsm, and gene expression in the .X
    
    Parameters
    ----------
    adata: anndata object
        10X provided anndata object
    drop_totalseq: bool
        If true, removes '_TotalSeq' from the end of .obsm[data_key] values
    data_key: str
        Name of data in .obsm to convert 
        
    Returns 
    -------
    adata: anndata object
        Dataframe with 'Antibody Capture' data moved to .obsm[data_key]
    '''
    adata.obsm[data_key] = pd.DataFrame(index=adata.obs_names)
    for i in adata.var_names[adata.var.feature_types=='Antibody Capture']:
        adata.obsm[data_key][i] = adata.obs_vector(i)
    adata = adata[:,adata.var.feature_types!='Antibody Capture'].copy()
    if drop_totalseq:
        adata.obsm[data_key].columns = pd.Series(adata.obsm[data_key].columns).str.split('_TotalSeq',expand=True)[0].values
    return adata

def default_kwargs(kwargs: dict, defaults: dict) -> dict:
    '''
    Helper function to combine lists of kwargs with defaults, if no kwargs is specified will use defaults instead
    '''
    for key in defaults.keys():
        if key not in kwargs.keys():
            kwargs[key] = defaults[key]
    return kwargs

def intersect_features(adatas: List[anndata.AnnData], data_key: str = 'protein') -> List[anndata.AnnData]:
    '''
    Subsets each anndata object to only genes and values of data_key found in every anndata objects
    
    Parameters
    ----------
    adatas: anndata object
        Each anndata object must have gene names in .var_names and must have the specified data_key in .obsm
    data_key: str
        Each anndata object's .obsm must have data_key as one of its fields
        
    Returns 
    -------
    adatas: anndata object
        List of anndata objects with same fields as inputted object, but containing information on only genes and 'data_key' that exists in every object
    '''
    intersection_GEX = list(set.intersection(*(set(val) for val in [adata.var_names for adata in adatas])))
    if not data_key is None:
        intersection_ADT = list(set.intersection(*(set(val) for val in [adata.obsm[data_key].columns for adata in adatas])))
    for i,adata in enumerate(adatas):
        adata = adata[:,intersection_GEX].copy()
        if not data_key is None:
            adata.obsm[data_key] = adata.obsm[data_key][intersection_ADT].copy()
        adatas[i] = adata
    return adatas

def preprocess_adatas(adatas: Union[anndata.AnnData, List[anndata.AnnData], str, List[str]] = None, 
                      convert_from_10X: bool = True,
                      make_unique: bool = True,
                      drop_totalseq: bool = True, 
                      intersect: bool = True, 
                      backup_urls: Union[str, List[str]] = None,
                      data_key: str = 'protein',
                      log_CP_GEX: int=1e4,
                      log_CP_ADT: int=1e2):
    
    '''
    Function to load and preprocess adatas from either filename or backup_url
    
    Parameters 
    ----------
    adatas: anndata or list of anndatas
        Individual or list of filepaths or anndata objects 
    convert_from_10X: bool
        Determines whether to convert from 10x data style to TODO
    make_unique: bool
        Determines whether to make labels in anndata objects unique by adding '-1' etc. to duplicates
    drop_totalseq: bool
        Determines whether to drop '_TotalSeq' from the data_key labels
    intersect: bool
        Determines whether to subset each dataframe to 'data_key' and gene expression data only found on all dataframes
    backup_urls: str or list of str
        A string or list of strings of urls to the anndata objects to be used if the filepath does not work
    data_key: str
        Must be in .obsm, name used for data intersection and adt normalization
    log_CP_GEX: bool or int
        Normalizes Gene expression data to log of 1 + counts per this value, default 1000
    log_CP_ADT: bool or int
        Normalizes ADT data to log of 1 + counts per this value, default 10
        
    Returns
    -------
    adatas
        List of anndata object(s)
            if make unique:
                containing unique obs and var labels, with duplicates named '-1', '-2' etc. 
            if convert_from_10x:
                with protein data moved to the .obsm
            if intersect:
                containing only the subset of genes and 'data_key' that exists in all objects
            if log_CP_GEX:
                log normalized gene counts per x value (default 1000)
            if log_CP_ADT: 
                log normalized 'data_key' counts per x value (default 10)
     '''
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
    adata: anndata object
    batch_key: str
        Column within the adata.obs that delineates batches
    sort: bool
        Whether to sort the batch names according to the sorted() function
    
    Returns 
    -------
    batch_masks: list of list of bool
        Boolean masks for each of the batches in batch_key
    batches: list of str
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

def obsm_to_X(adata: anndata.AnnData, data_key: str='protein') -> anndata.AnnData:
    '''
    Makes the .obsm[data_key] of an anndata object into its .X
    
    Parameters
    ----------
    adata: anndata object
    data_key: str
        key in the adata.obsm to convert to the .X of this new adata
    Returns 
    -------
    adata_new: anndata object
        Object with .X of anndata.obsm[data_key]
    '''
    adata_new = anndata.AnnData(adata.obsm[data_key],adata.obs,pd.DataFrame(index=adata.obsm[data_key].columns),dtype=float)
    return adata_new

def _marker(keyword: str, cols: List[str], allow_multiple: bool = False) -> Union[List[str], str]:
    '''
    Function to search for a marker in a list or listlike of markers
    
    Parameters
    ----------
    keyword: str
        Name of the marker you wish to search for
    cols: Listlike of str
        List of marker names in which to search for [keyword]
    allow_multiples: bool
        Whether to return a list of matches if there are multiple mathches or throw an error
        
    Returns 
    -------
    if allow_multiples 
        marker_names: List of str
            Names of the marker in cols
    else 
        marker_name: str
            Name of the marker in cols
    '''
    keyword = keyword.lower()
    marker_name = [i for i in cols if keyword in i.lower()]
    if marker_name == []: # If no marker was found
        raise ValueError("\"" + keyword + "\" not found in marker list.") 
    elif allow_multiple: # If we're allowing multiple just shove out whatever
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
        raise ValueError(f" '{keyword}' found multiple times in marker list: {marker_name} ") 
    else:
        return marker_name[0]

def marker(adata: anndata.AnnData, parameter: Union[str, List[str]], 
           data_key: Optional[str]='protein', allow_multiple: bool = False) -> Union[str, List[str]]:
    '''
    Lookup a marker name within the .X or the .obsm[data_key]
    
    Parameters
    ----------
    adata: anndata object
        Object with data in adata.X or adata.obsm[data_key]
    parameter: str or list of str
        Parameters to search for in adata
    data_key: optional, str
        Location in the .obsm to look for marker name
        If None, searches adata.var.index
    allow_multiple: bool
        Whether to return a list of matches if there are multiple mathches or throw an error
        
    Returns 
    -------
    if allow_multiples 
        marker_names: List of str
            Names of the marker in the adata object
    else 
        marker_name: str
            Name of the marker in the adata object
    '''
    
    if not data_key is None:
        cols = adata.obsm[data_key].columns.to_list()
    else:
        cols = adata.var.index.to_list()
    
    if isinstance(parameter, list): # Run recursively if it's a list of inputs
        results = []
        for i, word in enumerate(parameter):
            results.append(_marker(word, cols, allow_multiple))
        if allow_multiple:
            return sum(results,[])
        else:
            return results    
    else:
        return _marker(parameter, cols, allow_multiple)
        
def umap_thresh(adata: anndata.AnnData, h,
                markers: List[str]=None, batch_key: str='donor',
                data_key: Optional[str]='protein', umap_basis: str='X_umap',
                cmap: Optional[List[str]] =None ,**kwargs):
    '''
    Creates a umaps for the listed markers with expression data for the markers overlayed on top. 
    This visualization can be useful to see where protein and genes are expressed in the umap
    
    Parameters
    ----------
    adata: anndata object
    h: hierarchy object
    markers: list of str
        Markers to create umaps for
        If None (default), creates umap for all markers 
    batch_key: str
        Name of batch in adata to use for umap, uses batch_iterator to find key
    data_key: optional, str
        Location in the .obsm to look for marker name
        If None, searches adata.var.index
    umap_basis: str
        Passed to scanpy.pl.embedding as basis
    cmap: list of colors to use or matplotlib colormap
        Color map information used to color the gene expressions on the umap
    kwargs: dict
        Arguments to be passed to scanpy.pl.embedding
    
    '''
    if markers is None:
        all_markers = []
        for classification in h.get_classifications():
            for mm in h.tree[classification].data.markers:
                all_markers.append(mm)
        markers = list(set(all_markers))
    

    for mm in markers:
        if mm.endswith("_lo") or mm.endswith("_hi"):
            m = mm[:-3]
        else:
            m = mm
        to_del = False
        try:
            m = marker(adata,m,data_key)+"_mod_"+data_key
            adata.obs[m] = adata.obsm[data_key][m.split('_mod_')[0]]
            to_del = True
        except:
            m = marker(adata,m.split('_gex')[0],None)
        data = get_data(adata,m,data_key)
        for mask,b in zip(*batch_iterator(adata,batch_key)):
            t = h.get_threshold_info(mm,None,b,flexible_batch=True)[1]
            t = sorted(t)
            t[-1] = min(t[-1],max(data))
            t[0] = min(t[0],max(data))
            t[0] = max(t[0],min(data))
            t[-1] = max(t[-1],min(data))
            if not type(cmap) is str:
                from matplotlib import cm
                from matplotlib.colors import ListedColormap
                n_o = np.ceil(2560*(max(data)-max(t))/max(data))
                n_g = np.ceil(2560*(max(t)-min(t))/max(data))
                n_b = np.ceil(2560*min(t)/max(data))
                orange = cm.get_cmap('Oranges')
                blue = cm.get_cmap('Blues')
                green = cm.get_cmap('Greens')
                newcolors =  list(blue(np.linspace(0, 1, int(n_b)))[int(np.ceil(0.3*n_b)):int(np.ceil(0.860*n_b))]) + \
                             list(green(np.linspace(0, 1, int(n_g)))[int(np.ceil(0.3*n_g)):int(np.ceil(0.860*n_g))]) + \
                             list(orange(np.linspace(0, 1, int(n_o)))[int(np.ceil(0.3*n_o)):int(np.ceil(0.860*n_o))])
                newcolors = np.insert(newcolors,0,np.array([.95,.95,.95,1]),axis=0)
                cmap = ListedColormap(newcolors)

            sc.pl.embedding(adata[mask],basis = umap_basis,color=m,cmap=cmap,vmax=max(data),vmin=0,title=f'{b} {mm}',**kwargs)
        if to_del:
            del adata.obs[m]
    return
        
def get_data(adata: anndata.AnnData,
             parameter: str,
             preferred_data_key: Optional[str] = None,
             return_source: bool = False,
             force_preferred: bool = False) -> np.array:
    '''
    Searches an AnnData object along its .var, and .obs for a specified parameter
    To search the .obsm or .layers, define an appropriate data_key
                ***Adapted from farberanalyze package***
        
    Parameters
    ----------
    adata: AnnData object
        Object to search through     
    parameter: str
        string to search the AnnData object for. Append:
        "_gex" to limit checking to the .var,
        "_obs" to limit checking to the .obs,
        "_mod_" followed by an obsm key/layer name to force checking there.
        
        Note, these symbols overwrite preferred data key, and it is assumed
        that these symbols are not used elsewhere for variable names
    preferred_data_key: str
        Used to point to a key in .layers or in .obsm to check in
    return_source: bool
        Whether to return where the parameter was found in addition to the list of matches
    force_preferred: bool
        Whether to raise an error if the parameter is not found in the preferred_data_key
    
    Returns
    -------
    data: np array
       Array of data in the anndata object whoes label matches parameter
    '''
    data = None
    if "_mod_" in parameter:
        force_preferred = True
        preferred_data_key = parameter.split('_mod_')[1]
        parameter = parameter.split('_mod_')[0]
    if "_gex" in parameter or "_obs" in parameter:
        preferred_data_key = None

    if preferred_data_key is not None:
        assert preferred_data_key not in adata.obsm.keys() or preferred_data_key not in adata.layers.keys(), \
        f'Found {preferred_data_key} in both adata.obsm and adata.layers, please inspect AnnData Object and rectify.'
        
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
            return data, parameter+'_mod_'+preferred_data_key
        else:
            return data
    
    if force_preferred:
        raise ValueError(f'Did not find {preferred_data_key} in adata.obsm or adata.layers\n'+
                         f'Or did not find {parameter} in that data')
    else:
        assert parameter not in adata.var.index or parameter not in adata.obs.keys(), \
        'Found {parameter} in both adata.obs and adata.var, please inspect AnnData Object and rectify.'

        try: # Try checking .X
            try:
                if parameter.split("_gex")[0] in adata.var.index:
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

def list_tools(list1: list, operator: str,
               list2: list=[]) -> list:
    '''
    Helper functions to combine lists. Takes in two lists, with an operator for how to join the lists
    
    Operator can be
        '==' followed by a value - Returns a boolean list of whether each value is equal to the value after ==
        '+'  - Returns the sum of the lists
        '-'  - Returns the values in list1 that are not in list2
        '*'  - Returns the all permutations of sum of one element from list1 and one element from list2, in the order of 
               [list1[0] + list2[0], list1[1] + list2[0], list1[2] +list2[0], ... , list1[n] + list2[m]]
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
    
def generate_exlcusive_features(adata_list: Union[List[str], List[anndata.AnnData]], 
                                data_key: str = None):
    '''
    Reads in adata objects from a list of paths, or takes in a list of adata objects
    and finds the features in all objects
    
    Parameters
    ----------
    adata_list: listlike of anndata objs or strings
        List of objects in which to search for common features
    data_key: str
        Optional key in the adata.obsm to use for features, will be interesected in addition to .var_names
        
    Returns
    -------
    features: list of str
        features that appear in all anndata objects
    '''
    features = []
    for adata in adata_list:
        if isinstance(adata, str):
            adata = anndata.read_h5ad(adata_path)
        if data_key is None:
            features_one = list(adata.var_names)
        else:
            features_one = list(adata.var_names) + list(adata.obsm[data_key].columns)
        if features != []:
            features = [f for f in features_one if f in features]
            logg.print("intersection taken")
        else:
            logg.print("No itersection attempted")
            features = features_one
    return features