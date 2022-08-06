import anndata
import sys
import os
import pandas as pd
import numpy as np
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable
import scanpy as sc
from .logger import logg

def convert_10x(adata, drop_totalseq = True, data_key = 'protein'):
    '''TODO docstring, and add in logging.'''
    adata.obsm[data_key] = pd.DataFrame(index=adata.obs_names)
    for i in adata.var_names[adata.var.feature_types=='Antibody Capture']:
        adata.obsm[data_key][i] = adata.obs_vector(i)
    adata = adata[:,adata.var.feature_types!='Antibody Capture'].copy()
    if drop_totalseq:
        adata.obsm[data_key].columns = pd.Series(adata.obsm[data_key].columns).str.split('_TotalSeq',expand=True)[0].values
    return adata

def intersect_features(adatas, data_key='protein'):
    '''TODO docstring, and add in logging.'''
    intersection_GEX = list(set.intersection(*(set(val) for val in [adata.var_names for adata in adatas])))
    if not data_key is None:
        intersection_ADT = list(set.intersection(*(set(val) for val in [adata.obsm[data_key].columns for adata in adatas])))
    for i,adata in enumerate(adatas):
        adata = adata[:,intersection_GEX].copy()
        if not data_key is None:
            adata.obsm[data_key] = adata.obsm[data_key][intersection_ADT].copy()
        adatas[i] = adata
    return adatas

def preprocess_adatas(adatas=None, convert_from_10X = True, make_unique = True, drop_totalseq = True, 
                      intersect=True, backup_urls = None, data_key = 'protein',log_CP_GEX=True,log_CP_ADT=True):
    ''' Function to load and preprocess adatas from either filename/backup_url strings,  TODO docstring, and add in logging.'''
    if isinstance(adatas, anndata.AnnData) or isinstance(adatas,str):
        adatas = [adatas]
        backup_urls = [backup_urls]
    if backup_urls is None:
        backup_urls = [None] * len(adatas)
    for i,(adata,backup_url) in enumerate(zip(adatas,backup_urls)):
        if isinstance(adata,str):
            if backup_url is None:
                adata = anndata.read_h5ad(adata)
            else:
                adata = sc.read_10x_h5(adata, backup_url=backup_url, gex_only=False)
                
        if make_unique:
            adata.var_names_make_unique()
            adata.obs_names_make_unique()
            
        if convert_from_10X:
            adata = convert_10x(adata,drop_totalseq)
        
        if drop_totalseq and make_unique and convert_from_10X:
            assert len(adata.obsm[data_key].columns) == len(adata.obsm[data_key].columns.unique()), 'Removing TotalSeq tags made antibodies not unique, try again without `drop_totalseq = True`'
        adatas[i] = adata
        
    if intersect:
        adatas = intersect_features(adatas,data_key)
        adatas = adatas
    for i, adata in enumerate(adatas):
        if log_CP_GEX > 0 and (adata.X.astype(int) != adata.X).nnz==0:
            if not isinstance(log_CP_GEX, int):
                log_CP_GEX = 1e4
            adata.layers['counts'] = adata.X.astype(int).copy()
            sc.pp.normalize_total(adata, target_sum = log_CP_GEX)
            sc.pp.log1p(adata)
        if (log_CP_ADT > 0) and (data_key in adata.obsm.keys()) and all(adata.obsm[data_key].astype(int) == adata.obsm[data_key]):
            if not isinstance(log_CP_ADT, int):
                log_CP_ADT = 1e2
            adata.obsm[f'{data_key}_counts'] = adata.obsm[data_key].astype('int32').copy()
            adata.obsm[data_key] = (adata.obsm[data_key].div(adata.obsm[data_key].sum(axis=1)/log_CP_ADT, axis=0))
            adata.obsm[data_key].apply(np.log1p).astype('float32')
            adata = adata[adata.obsm[data_key].isna().sum(axis=1)==0]
        adatas[i] = adata
    return adatas


def batch_iterator(adata, batch_key,sort=True):
    '''
    adata: anndata object
    batch_key: str, column within the adata.obs that delineates batches
    returns a list of boolean masks and a corresponding list of batch names that lines up to it
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

def obsm_to_X(adata, data_key='protein'):
    ''' Creates a new anndata object by copying the anndata
    adata: anndata object
    data_key: key in the adata.obsm to convert to the .X of this new adata
    returns new anndata object
    '''
    adata_new = anndata.AnnData(adata.obsm[data_key],adata.obs,pd.DataFrame(index=adata.obsm[data_key].columns))
    return adata_new


def _marker(keyword, cols, allow_multiple = False):
    ''' Internal function for marker name lookup '''
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

def marker(adata, parameter, data_key = 'protein', allow_multiple = False):
    ''' Lookup a marker name within the .X or the .obsm[data_key]
    adata: anndata object with data in adata.X or adata.obsm[data_key]
    parameter: str or list of str of parameters to look up
    data_key: str, if not None, location in the .obsm to look for marker name.
    allow_multiple: bool, whether to return a list of matches, or throw an error.
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
        

def get_data(adata: anndata,
             parameter: str,
             preferred_data_key: Optional[str] = None,
             return_source = False,
             force_preferred = False):
    '''
    - adapted from farberanalyze package
    Searches an AnnData object along its .var, and .obs for a specified parameter
    To search the .obsm or .layers, define an appropriate data_key
        
    Parameters
    ----------
    adata
        AnnData     
    parameter
        string to search the AnnData object for. Append:
        "_gex" to limit checking to the .var,
        "_obs" to limit checking to the .obs,
        "__" followed by an obsm key/layer name to force checking there.
        Note, these symbols overwrite preferred data key, and it is assumed
        that these symbols are not used elsewhere for variable names
    preferred_data_key
        used to point to a key in .layers or in .obsm to check in
        
    Returns
    -------
    A list of values
    '''
    data = None
    
    if "__" in parameter:
        force_preferred = True
        preferred_data_key = parameter.split('__')[1]
        parameter = parameter.split('__')[0]
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
            return data, parameter+'__'+preferred_data_key
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
                data = np.array(adata.obs[parameter.split("_obs")])
            except:
                parameter = _marker(parameter.split("_obs"),adata.obs.columns.to_list(), False)
                data = np.array(adata.obs[parameter])
        except:
            pass
        
        if not data is None:
            if return_source:
                return data, parameter.split("_obs")+'_obs'
            else:
                return data

    if preferred_data_key is None:
        raise ValueError(f'Did not find {parameter} in adata.var or adata.obs') 
    else:
        raise ValueError(f'Could not find {parameter} in adata.')
    return

def list_tools(list1,operator,list2=[]):
    '''Helper functions to combine lists. Takes in two lists, with an operator for how to join the lists'''
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
    
def generate_exlcusive_features(adata_list, data_key = None):
    ''' Reads in adata objects from a list of paths, or takes in a list of adata objects.
    Returns the intersection of possible features from each adata in the list. 
    adata_list: iterable of anndata objects or str
    data_key: optional key in the adata.obsm to use for features'''
    features = []
    for adata in adata_list:
        if isinstance(adata, str):
            adata = anndata.read_h5ad(adata_path)
        else:
            if data_key is None:
                features_one = list(adata.var_names)
            else:
                features_one = list(adata.var_names) + list(adata.obsm[data_key].columns)
            if features != []:
                print("intersection taken")
                features = [f for f in features_one if f in features]
            else:
                print("No itersection attempted")
                features = features_one
    return features