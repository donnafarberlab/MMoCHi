import numpy as np
import pandas as pd
import sklearn 
import sklearn.ensemble
import sklearn.neighbors
from sklearn.calibration import CalibratedClassifierCV

import imblearn
import random
import gc
from collections import Counter
import scipy.sparse as sp
import scanpy as sc
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable, List, Set, Iterable

from . import hierarchy
from . import utils
from .logger import logg
from . import thresholding


DEBUG_ERRORS = False

def classifier_setup(adata: anndata.AnnData, x_modalities: Union[str,List[str]],
                     data_key: Optional[str]=utils.DATA_KEY, reduce_features_min_cells: int=0, 
                     features_limit=None) -> Tuple[sp.csr_matrix, list]:
    """
    Setup that can optionally be completed before running mmc.classify. This can be run before the classifier (to reduce runtime of the classifier function in a parameter optimization loop) or is automatically run when training a classifier. It concatenates the .X and any data_key in the .obsm, then performs feature reduction (if reduce_features_min_cells > 0). Next, features can be limited by an external feature set. Then, it sorts the resulting feature_names (the columns from the .X and .obsm[data_key]) and csr.matrix, alphabetically, to make the feature order reproducible across runs. If defined, feature limits can be performed so that you can match the expected features of the hierarchy.

    Parameters
    ----------
    adata
        Object containing gene expression data, and expression data for modalities for every data key in `.obsm`
    x_modalities
        Name of the modality of the data in the .X of adata
    data_key
        Key in adata.obsm to concatenate into .X and to reduce features across
    reduce_features_min_cells
        Remove features that vary in fewer than this number of cells passed to _reduce_features
    features_limit
        listlike of str or dictionary in the format `{'modality_1':['gene_1', 'gene_2', ...], 'modality_2':'All'}`
        Specifies the allowed features to classify on for a given modality
    
    Returns
    -------
    scipy.sparse.csr_matrix
        Reduced adata data for classification
    list
        List of features that were checked/used in the reduction process 
    """
    logg.print('Setting up...')
    if type(x_modalities) is str and x_modalities in adata.var.columns:
        x_modalities = adata.var[x_modalities]
    if data_key is None:
        logg.print('Using .X only')
        X = sp.csr_matrix(adata.X)
        utils._validate_names(adata.var_names) 
        features_used = list(adata.var_names+"_mod_"+x_modalities)
    else:
        logg.print('Using .X and '+str(data_key))
        X = sp.csr_matrix(adata.X)
        X = sp.hstack([X,sp.csr_matrix(adata.obsm[data_key])],format='csr')
        utils._validate_names(adata.var_names) 
        utils._validate_names(adata.obsm[data_key].columns) 
        features_used = list(adata.var_names+"_mod_"+x_modalities) + list(adata.obsm[data_key].columns+"_mod_"+data_key)
    
    
    assert len(features_used) == len(set(features_used)), 'All feature names must be unique...'
    
    if type(features_limit) is str and features_limit in adata.var.columns:
        logg.warn(f"Limiting .X features using {features_limit} column in adata.var...")
        features_limit = adata.var_names[adata.var[features_limit]] +"_mod_"+ x_modalities
        features_limit = _limit_features_list_to_dict(features_limit)
        if not data_key is None:
            features_limit[data_key] = 'All'
            
    X, features_used = _reduce_features(X,features_used, min_cells=reduce_features_min_cells)
    X, features_used = _limit_and_realign_features(X, features_used, features_limit)
    return X,list(features_used)

def _limit_features_list_to_dict(features_limit: List[str]) -> dict:
    """
    Creates dictionary of feature names to modalities, splits rows of np.array with first column as key and second column as value.
    
    Parameters
    ----------
    features_limit
        array with feature names in first column and modalities in second column
    
    Returns
    -------
    A dictionary with feature names as keys and modalities as values
    """
    features_limit = np.array(features_limit)
    features_limit_mods = np.array([i.split('_mod_')[-1] for i in features_limit])
    features_limit_names = np.array([i.split('_mod_')[0] for i in features_limit])
    f_limit = dict()
    for i in set(features_limit_mods):
        f_limit[i] = features_limit_names[features_limit_mods==i]
    return f_limit

def _limit_and_realign_features(X: sp.csr_matrix, feature_names: List[str],
                                features_limit: Union[List[str], dict],
                                assert_all_included: bool=False) -> Tuple[sp.csr_matrix, list]:
    """
    Removes any features from X that are not in features_limit (if supplied) and realigns new matrix of features and features list to reflect removed features, in addition to alphabetically sorting the features.
    
    Parameters
    ----------
    X
        sparse matrix containing feature information
    feature_names
        names of features in X
    features_limit
        features to keep in X
    assert_all_included
        Whether to throw error if not all features in features limit also exist in features names
    Returns
    -------
    check_array: csr.spare_matrix
        Validated, limited, realigned, and sorted array
    features_used: list
        The sorted features in check_array
    """
    if (features_limit is None) or (len(feature_names) == len(features_limit)) & (list(feature_names) == list(features_limit)):
        features_used = feature_names
    else:
        feature_names = np.array(feature_names)
        if type(features_limit) is str:
            features_limit = {features_limit:'All'}
        if type(features_limit) is dict:
            logg.warn('Converting dict of features_limit to features by modality')
            f_limit = list()
            for i in features_limit:
                if type(features_limit[i]) is str and features_limit[i] == 'All':
                    new_adds = feature_names[np.char.endswith(feature_names,'_mod_'+i)]
                else:
                    new_adds = [j+ '_mod_' + i for j in features_limit[i]]
                f_limit.extend(new_adds)
        else:
            logg.warn('Limiting features using listlike, make sure you include modality information...')
            f_limit = features_limit
        assert pd.api.types.is_list_like(f_limit), f'features_limit of type {type(f_limit)} not supported'
        feature_names = np.array(feature_names)
        feature_mask = np.isin(feature_names, f_limit)
        logg.debug(f'Limiting features, from {len(feature_names)} to {sum(feature_mask)}')
        assert sum(feature_mask) > 0, 'Limiting to 0 features invalid, check input...'
        X=X[:,np.flatnonzero(feature_mask)]
        features_used = feature_names[np.flatnonzero(feature_mask)]
    if sorted(features_used) != list(features_used):
        logg.info('Resorting to enforce sorted order of features by name')
        sort = np.argsort(features_used)
        X = X[:,sort]
        features_used = np.array(features_used)[sort]
        logg.debug(f'Features_used are:{features_used}')
    if assert_all_included:
        assert len(features_used) == len(features_limit), \
        f'Not all features in feature_limit availible in feature_names {set(features_limit) - set(features_used)}'
    return sklearn.utils.check_array(X,accept_sparse=True), list(features_used)    

def _reduce_features(X: sp.csr_matrix, feature_names: List[str],
                     min_cells: int=50) -> Tuple[sp.csr_matrix, List[str]]:
    """
    Remove features that are expressed in fewer than min_cells_threshold cells.
    
    Parameters
    ----------
    X
        Contains feature information
    feature_names
        Contains feature names
    min_cells
        Minimum cells threshold. Any feature expressed in fewer cells will be excluded.
        
    Returns
    -------
    X: 2D scipy.sparse.csr_matrix
        Contains feature information for features expressed in at least min_cells
    feature_names: listlike 
        Contains feature names for features expressed in at least min_cells
    """
    if min_cells > 0:
        num_cells = X.getnnz(axis=0)
        logg.info(f"Removing {sum(num_cells<min_cells)} features lacking expression in a minimum of {min_cells} events...")
        X=X[:,np.flatnonzero(num_cells>=min_cells)]
        feature_names = np.array(feature_names)[np.flatnonzero(num_cells>=min_cells)]
    else:
        logg.info("No min_cells feature reduction...")
    return X, feature_names


def classify(adata: anndata.AnnData, hierarchy: hierarchy.Hierarchy, 
             key_added: str='lin', data_key: Optional[str]=utils.DATA_KEY,
             x_data_key: Optional[str]=None, x_modalities: str='GEX',
             batch_key: Optional[str]=utils.BATCH_KEY, retrain: bool=False,
             plot_thresholds: bool=False, reduce_features_min_cells: int=25,
             allow_spike_ins: bool=True, enforce_holdout: bool=True, 
             probabilities_suffix: Optional[str]='_probabilities', 
             resample_method: Optional[str]='oversample', weight_integration: bool=False,
             X: Optional[sp.csr_matrix]=None, features_used: Optional[List[str]]=None,
             probability_cutoff: float=0,
             features_limit: Optional[Union[str, List[str], dict]]=None, 
             skip_to: Optional[str]=None, end_at: Optional[str]=None,
             reference: Optional[str]=None) -> Tuple[anndata.AnnData, hierarchy.Hierarchy]:
    """
    Classify subsets using the provided hierarchy. If the hierarchy has not yet been trained, this function trains the hierarchy and predicts cell types. If the hierarchy has been trained, this classifier will not retrain random forests unless explicitly asked to (retrain=True). Many classifier kwargs are defined in the hierarchy on a node-by-node basis. Please see mmc.Hierarchy() and mmc.add_classification() for additional details. 
    
    Classifications information for every level are recorded in .obsm['lin']. [classification_layer]_class denotes the classified MMoCHi classification for each event. [classification_layer]_hc denotes the high-confidence thresholded categories for each event. [classification_layer]_hold_out denotes whether the event is part of explicit hold-out for that layer for each event. [classification_layer]_train is a boolean column indicating whether the event was used in training (noise events can be neither hold_out nor _train). [classification_layer]_traincounts denotes the number of times an event was used in training (amplified events can have _traincounts >1). If SMOTE or its variations are used to balnce training data, traincounts will refer to rebalanced training counts and will not account for generated training events. [classification_layer]_probability indicates the confidence of the classifier at the specified layer for each event.
    
    Note: Only [classification_layer]_hc and [classification_layer]_class layers are created for cutoff layers
    
    We recommend following this function with mmc.terminal_names() to put classifications into a `.obs` column
    
    Note: If this function terminates prematurely, adata.obs_names may be corrupted. The original adata.obs_names can be found in adata.obs['MMoCHi_obs_names'].

    Parameters
    ----------
    adata
        Object containing expression of one modality in the .X and optionally expression data for another modality in the .obsm[data_key]
    hierarchy
        Hierarchy design specifying one or more classification levels and subsets.
    key_added
        Key in adata.obsm to store the resulting DataFrame of the classifier's results
    data_key
        Key in adata.obsm to be used for high-confidence thresholding.
    x_data_key
        Key in adata.obsm used for predictions, if undefined, defaults to data_key
    x_modalities
        Column in adata.var (if x_modalities is in adata.var.columns), or name of the modality of the data in the .X
    batch_key
        Name of a column in adata.obs that corresponds to a batch for use in the classifier   
    retrain
        If the classification level of a hierarchy is already defined, whether to replace that classifier with a retrained one.
        If False, no new classifier would be trained, and that classifier is used for prediction of the data.
        Overrides reduce_features_min_cells to 0 (as this early-stage feature reduction often breaks re-running models).
    plot_thresholds
        Whether to display histograms of expression distribution while performing high-confidence thresholding
    reduce_features_min_cells
        Remove features that are expressed in fewer than this number of cells, passed to _reduce_features. 
        Feature reduction can be a very powerful tool for improving classifier performance. 
        Only used if X and features_used are not provided and retrain is True.
    allow_spike_ins
        Whether to allow spike ins when doing batched data. Spike ins are performed if a subset is below the minimum threshold within an individual batch, 
        but not the overall dataset. If False, errors may occur if there are no cells of a particular subset within a batch. Warning, spike ins currently 
        do not respect held out data, meaning there will be much fewer than 20% held out data if spike ins are frequent.
    enforce_holdout
        Whether to enforce a hold out of an additional 20% at each classification level, to prevent spike-ins from using that data.     
    probabilities_suffix
        If defined, probability outputs of the full df of predict_probabilities will be saved to the .uns[level+probabilities_suffix]
    resample_method
        Method used to resample high-confidence data prior to training, passed to _balance_training_classes. "Oversample" or None. Other methods are experimental. See _balance_training_classes and the imblearn package for details.
    weight_integration
        Whether to have each batch represented equally in the forest, or weight the number of trees from each batch to their representation in the total dataset 
    X, features_used scipy.sparse.csr_matrix and listlike of str
        Optional setup data to circumvent the initial classifier setup functions. The presence of both of these overrides any predefined feature reduction 
        techniques.   
    probability_cutoff
        Between 0 and 1, the minimum proportion of trees that must agree on a call for that event to be labeled by the classifier. This is experimental.
    features_limit
        A list of str specifying features to limit to, in the format [feature_name]_mod_[modality], e.g. 
    skip_to
        Name of hierarchy classification levels to start at. This is extremely useful when debugging particular levels of the classifier.
    end_at
        Name of hierarchy classification levels to prematurely end at. This is extremely useful when debugging particular levels of the classifier.
        The order of levels corresponds to the order they are defined in the hierarchy. Requires that the AnnData object has a predefined .obsm[key_added]. If
        skip_to is defined, it replaces only the columns that occur at and beyond that level.
    reference
        Column in the adata.obs to be used for comparing (in the log file) the results of high-confidence thresholding to predetermined annotations or clusters
        
    Returns
    -------
    adata: AnnData object
        Object containing a .obsm[key_added] containing columns corresponding to the classification high-confidence ("_hc"), the held out data not used for testing ("_hold_out"), the prediction ("_class") for each level. If probababilities_suffix is defined, also contains "_probabilities" keys in the .uns corresponding to the prediction probabilities of each class for each event.
    Hierarchy: Hierarchy object
        A trained hierarchy with classifiers built into it and a record of thresholds and other settings used. 
    """
    if x_data_key is None:
        x_data_key = data_key
    setup_complete = (not X is None) and (not features_used is None)
    if not setup_complete:
        if retrain is False:
            reduce_features_min_cells = 0
        X,features_used = classifier_setup(adata,x_modalities,x_data_key,reduce_features_min_cells,features_limit=features_limit)
    features_used_X = features_used
    logg.print(f"Set up complete.\nUsing {len(features_used_X)} features")
    classifications = hierarchy.get_classifications()
    if skip_to is None:
        adata.obsm[key_added] = pd.DataFrame(index=adata.obs.index)
        adata.obsm[key_added]['All_class'] = 'All'
    else:
        assert skip_to in classifications, 'Attempting to skip to an invalid classification level'
        classifications = classifications[classifications.index(skip_to):]
        parent,gparent = hierarchy.classification_parents(skip_to)
        assert gparent+"_class" in adata.obsm[key_added].columns, 'Attempting to skip to a level beneath an unclassified level'
        logg.warning(f'Skipped to classification on {skip_to}')
        for classification in classifications:
            adata.obsm[key_added].drop([classification+"_class",classification+"_hc",classification+"_hold_out",
                                        classification+"_probability",classification+"_train",classification+"_traincounts"],
                                        axis=1, inplace=True, errors='ignore')
    if not end_at is None:
        classifications = classifications[:classifications.index(end_at)+1]     
    total_adata, total_X = adata, X  
    total_adata.obs['MMoCHi_obs_names'] = total_adata.obs_names
    total_adata.obs_names = range(0,len(total_adata))    
    total_adata.obs_names = total_adata.obs_names.astype(str)
    
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    if weight_integration:
        weights = list(np.array([sum(b)/len(b) for b in batch_masks]))
    else:
        weights =  [1/len(batches)] * len(batches)
    logg.info(f'Using weights of: {weights} for random forest n_estimators')
    for level in classifications:
        try:
            parent,gparent = hierarchy.classification_parents(level)
            assert gparent+"_class" in total_adata.obsm[key_added].columns, f"Cannot subset on {parent} to train {level} as {gparent} has not been classified"
            subset_mask = total_adata.obsm[key_added][gparent+'_class'] == parent
            assert sum(subset_mask) > 0, f"No events were classified as {parent} in {gparent}_class"
            subset_adata, subset_X = total_adata[subset_mask],total_X[subset_mask.values] 
            features_used = features_used_X
            logg.print("Data subsetted on "+ parent+' in '+ gparent) 
            
            if retrain or not hierarchy.has_clf(level):
                logg.print(f'Running high-confidence populations for {level}...')         
                hc_dfs, hc_batches = [], []
                for batch_mask, batch in zip(batch_masks, batches):
                    adata = subset_adata[batch_mask[subset_mask]]
                    if not batch_key is None:
                        logg.info(f'Running high-confidence thresholds in {batch}')
                    if len(adata) > 0:
                        df = hc_threshold(adata,hierarchy,level,data_key,reference=reference,
                                          batch=batch,plot=plot_thresholds)
                        hc_dfs.append(df)
                        hc_batches.append(batch)
                        logg.debug(f"high-confidence thresholding on {len(df)} events")
                    else:
                        logg.info(f'No events in subsetted adata, skipping high-confidence thresholding')
                    del adata
                    gc.collect()

                if hierarchy.get_info(level,'is_cutoff'):
                    logg.print(f'Performing cutoff for {level}...')
                    df = pd.concat(hc_dfs)
                    df[level+'_class'] = df[level+'_hc']
                    df.loc[df[level+'_class'] == '?',level+'_class'] = np.nan
                else:
                    logg.print(f'Preparing training data for {level}...')
                    hc_df = pd.concat(hc_dfs)
                    min_events = hierarchy.get_info(level,'min_events')
                    
                    if min_events != 0:
                        logg.info(f'Checking subsets for minimum events...')
                        if min_events < 1:
                            min_events = np.ceil(min_events * len(subset_adata))
                            
                        subsets = hierarchy.subsets_info(level).keys()
                        vals = hc_df[level+'_hc'].value_counts()
                        vals.drop('?',errors='ignore',inplace=True)
                        logg.debug(f'Values before min_events ({min_events}) thresholding\n{vals.to_string()}')
                        logg.debug(f'len(vals) = {len(vals)}, len(subsets) = {len(subsets)}')
                        if min(vals) < min_events or len(vals) < len(subsets):
                            logg.debug(f'Must remove at least one subset')
                            removed_subsets = list(set(subsets)-set(vals.index[vals>=min_events]))
                            logg.warning(f'Removing {removed_subsets} cells as they are below {min_events} events')
                            hc_df.loc[hc_df[level+'_hc'].isin(removed_subsets),level+'_hc'] = '?'
                            subsets = list(set(subsets)-set(removed_subsets))
                        assert len(hc_df[hc_df[level+'_hc']!=0][level+'_hc'].unique()) > 1, 'Cannot train classifier with only one class'
                        
                        force_spike_ins = hierarchy.get_info(level,'force_spike_ins')
                        if enforce_holdout:# or hierarchy.get_info(level,'calibrate'):
                            original_hc_df = hc_df.copy()
                            logg.debug(f" original_hc_df {original_hc_df}")
                            hc_dfs_holdout = []
                            for subset in subsets:
                                hc_dfs_holdout.append(hc_df[hc_df[level+"_hc"]==subset].sample(frac=.9))
                            hc_df = pd.concat(hc_dfs_holdout)
                        if not batch_key is None:
                            hc_barcodes_batches = []
                            spike_in_barcodes = [] # saving spike_in barcodes to make sure they are all not hold out
                            for df, batch in zip(hc_dfs, hc_batches):
                                hc_barcodes_batch = list(df.index)
                                for subset in subsets:
                                    vals = df[level+'_hc'].value_counts()
                                    if not subset in vals.index or (vals.loc[subset] < min_events):
                                        if not subset in vals.index:
                                            events_needed = min(min_events, sum((hc_df[level+'_hc']==subset) & ~hc_df.index.isin(hc_barcodes_batch))-1) 
                                        else:
                                            events_needed = min(min_events - vals.loc[subset],sum((hc_df[level+'_hc']==subset) & \
                                                                                                  ~hc_df.index.isin(hc_barcodes_batch))-1) 
                                            events_needed = max(events_needed,0)
                                        if allow_spike_ins:
                                            logg.info(f'Spiking in {events_needed} of {subset} in {batch} to reach {min_events} events')
                                            spike_in_barcode_batch = random.sample(list(hc_df[(hc_df[level+'_hc']==subset) & \
                                                                                          ~hc_df.index.isin(hc_barcodes_batch)].index), events_needed)
                                            spike_in_barcodes += spike_in_barcode_batch
                                            hc_barcodes_batch += spike_in_barcode_batch
                                            
                                        else:
                                            assert False, f"{batch_key} of {batch} is missing subset {subset}. Cannot train without spike ins."
                                    if subset in force_spike_ins:
                                        events_needed = min(vals.loc[subset]*len(batches),sum((hc_df[level+'_hc']==subset) & \
                                                                                              ~hc_df.index.isin(hc_barcodes_batch))-1)
                                        events_needed = max(events_needed,0)
                                        logg.info(f'Spiking in {events_needed} of {subset} in {batch} as required by force_spike_ins)')
                                        spike_in_barcode_batch = random.sample(list(hc_df[(hc_df[level+'_hc']==subset) & \
                                                                                      ~hc_df.index.isin(hc_barcodes_batch)].index), events_needed)
                                        spike_in_barcodes += spike_in_barcode_batch
                                        hc_barcodes_batch += spike_in_barcode_batch
                                hc_barcodes_batches.append(hc_barcodes_batch)
                        else:
                            hc_barcodes_batches = [list(hc_df[hc_df[level+'_hc']!='?'].index)]
                    if enforce_holdout: # or hierarchy.get_info(level,'calibrate'):
                        hc_df = original_hc_df
                        
                    in_danger_noise_checker, max_training, clf_kwargs, class_weight, f_limit = \
                    hierarchy.get_info(level,['in_danger_noise_checker', 
                                              'max_training','clf_kwargs','class_weight', 'features_limit'])
                    logg.info(clf_kwargs)
                    if class_weight == 'balanced' and not batch_key is None:
                        y_org = hc_df.loc[hc_df[level+'_hc']!='?',level+'_hc'].values
                        # y_org = hc_df.loc[[True] * len(hc_df),level+'_hc'].values
                        classes = np.unique(y_org)
                        class_weight = sklearn.utils.class_weight.compute_class_weight('balanced', classes=np.unique(y_org), y=y_org)
                        class_weight = dict(zip(classes,class_weight))
                        class_weight = {}
                        for subset in classes:                            
                            # class_weight[subset]= math.log((1/(sum(y_org==subset)/len(y_org))+1),100)
                            class_weight[subset] = (1/(sum(y_org==subset)/len(y_org))) ** (1./2)
                        logg.info(f'Manually made a balanced class_weight: {class_weight}')
                    

                    logg.print(f'Initializing classifier for {level}...')
                    clf = sklearn.ensemble.RandomForestClassifier(class_weight=class_weight,
                                                                  warm_start=(not batch_key is None),**clf_kwargs)
                    subset_X, features_used = _limit_and_realign_features(subset_X, features_used, f_limit)
                    n_estimators = clf.n_estimators
                    clf.n_estimators = 0
                    dfs = []
                    train_dfs = []
                    for hc_barcodes_batch, batch, weight in zip(hc_barcodes_batches, hc_batches,weights):
                        clf.n_estimators += int(np.ceil(n_estimators * weight))
                        logg.info(f'Training {int(np.ceil(n_estimators * weight))} new estimators using {batch}...')
                        df = hc_df.loc[hc_barcodes_batch].copy()
                        mapper = dict(zip(subset_adata.obs_names,range(0,len(subset_adata.obs_names))))
                        hc_barcodes_batch = [mapper[i] for i in hc_barcodes_batch]
                        X,y, traincounts = _get_train_data(subset_X[hc_barcodes_batch],df,level,resample_method,max_training,in_danger_noise_checker, features_used)
                        train_dfs.append(pd.DataFrame.from_dict(dict(zip(df.index,traincounts)), orient="index", columns=["training"]))
                        clf.fit(X, y)
                        dfs.append(df)
                        del X, y
                    hierarchy.set_clf(level,clf,features_used)
                    df = pd.concat(dfs)
                    
                    lost_items = hc_df[~hc_df.index.isin(df.index)].copy()
                    lost_items[level+'_hold_out'] = True
                    lost_items[level+'_traincounts'] = 0
                    lost_items[level+'_train'] = False
                    df = pd.concat([df,lost_items])
                    
                    df.reset_index(inplace=True)
                    
                    df = df.groupby(['index',level+'_hc']).sum()
                    df.reset_index(level=1,inplace=True)
                    
                    train_dfs.append(pd.DataFrame(0,index=lost_items.index,columns=[level+'_training']))
                    train_df = pd.concat(train_dfs)
                    train_df.reset_index(inplace=True)
                    train_df = train_df.groupby(['index']).sum()
                    df[level+'_traincounts'] = train_df['training'].astype(int)
                    df[level+'_train'] = train_df['training'].astype(bool)
                    del train_df                
                    df[level+'_hold_out'] = (df[level+'_hold_out']) #& (df[level+'_hc'] != '?')

                logg.print(f"Merging data into adata.obsm['{key_added}']")
                
                total_adata.obsm[key_added] = total_adata.obsm[key_added].merge(df,how='left',left_index=True,
                                          right_index=True,suffixes=['_error',None]) 
                
                #spike ins might be holdout in one batch and trained with in another
                if not hierarchy.get_info(level,'is_cutoff') and not batch_key is None and len(spike_in_barcodes) > 0:
                    total_adata.obsm[key_added].loc[list(set(spike_in_barcodes)),level+'_hold_out'] = 0.0
                
                if not hierarchy.get_info(level,'is_cutoff'):
                    if hierarchy.get_info(level,'calibrate'):
                        logg.print(f'Running calibration on random forest')
                        
                        train_mask = total_adata[subset_mask].obsm[key_added][level+'_train'].astype(bool).values
                        hc_mask = (total_adata[subset_mask].obsm[key_added][level+'_hc'] != '?').values
                        third_mask = np.random.choice([True,False,False],len(hc_mask))
                        calibration_mask = ~train_mask & hc_mask & third_mask
                        # calibration_mask = train_mask
                        # if calibration_mask doesn't include all
                        y_calibration = total_adata[subset_mask].obsm[key_added][level+'_hc'][calibration_mask]
                        if ((sum(calibration_mask) < 1000) or 
                           (not set(y_calibration.values) == (set(total_adata[subset_mask].obsm[key_added][level+'_hc'].values) - {'?'}))):
                            calibration_mask = ~train_mask & hc_mask
                            y_calibration = total_adata[subset_mask].obsm[key_added][level+'_hc'][calibration_mask]
                            logg.info(f'Calibration will not include the 1/3rd mask')
                        y_calibration = y_calibration.values
                        # if len(y_calibration)>1000:
                        #     cal_method = 'isotonic'
                        # else:
                        #     cal_method = 'sigmoid'
                        cal_method='isotonic'
                        logg.info(f'Calibrating with method {cal_method}')
                        cal_clf = CalibratedClassifierCV(clf, method=cal_method, cv="prefit") # isotonic should perform better with imbalanced classes.
                        X_calibration = subset_X[calibration_mask]
                        
                        cal_clf.fit(X_calibration, y_calibration)
                        hierarchy.set_clf(level,cal_clf,features_used)
                        
            if not hierarchy.get_info(level,'is_cutoff'):
                df = pd.DataFrame(index=subset_adata.obs_names)

                logg.print(f'Predicting for {level}...')

                clf,f_names = hierarchy.get_clf(level)
                subset_X, features_used = _limit_and_realign_features(subset_X, features_used, f_names, True)
                assert list(f_names) == list(features_used), 'Order of f_names and features_used not matching. Aborting...'
                full_probabilities, df[level+'_probability'], df[level+'_class'] = _predict_probabilities_cutoff(subset_X, clf, probability_cutoff)

                if not probabilities_suffix is None:
                    total_adata.uns[level+probabilities_suffix] = pd.DataFrame(full_probabilities,
                                                                       index=np.array(total_adata.obs['MMoCHi_obs_names'])[df.index.astype(int)],columns=clf.classes_)

                logg.print(f"Merging data into adata.obsm['{key_added}']")
                total_adata.obsm[key_added] = total_adata.obsm[key_added].merge(df,how='left',left_index=True,
                                          right_index=True,suffixes=['_error',None])        
            logg.print('Predicted:')
            logg.print(df[level+'_class'].value_counts())
            if not reference is None:
                logg.debug(f"\n{pd.crosstab(df[level+'_class'],subset_adata.obs[reference]).T.to_string()}")

            del subset_X
            del subset_adata
            gc.collect()
        except Exception as e:
            logg.error(f'Ran into an error, skipping {level}',exc_info=1)
            if DEBUG_ERRORS:
                assert False
        # Attempt to convert to categories and bools to simplify
        
    del total_X
    cols = total_adata.obsm[key_added].columns
    try:
        logg.info(f'Converting columns in adata.obsm["{key_added}"] to savable dtypes...')
        total_adata.obsm[key_added][cols[cols.str.endswith('_class')]] = \
        total_adata.obsm[key_added][cols[cols.str.endswith('_class')]].astype(str).astype('category')
        total_adata.obsm[key_added][cols[cols.str.endswith('_hc')]] = \
        total_adata.obsm[key_added][cols[cols.str.endswith('_hc')]].astype(str).astype('category')

        total_adata.obsm[key_added].loc[:,cols.str.endswith('_hold_out')] = \
        total_adata.obsm[key_added].loc[:,cols.str.endswith('_hold_out')].fillna(False)
        total_adata.obsm[key_added][cols[cols.str.endswith('_hold_out')]] = \
        total_adata.obsm[key_added][cols[cols.str.endswith('_hold_out')]].astype(bool)

        total_adata.obsm[key_added].loc[:,cols.str.endswith('_train')] = \
        total_adata.obsm[key_added].loc[:,cols.str.endswith('_train')].fillna(False)
        total_adata.obsm[key_added][cols[cols.str.endswith('_train')]] = \
        total_adata.obsm[key_added][cols[cols.str.endswith('_train')]].astype(bool)
        total_adata.obsm[key_added].loc[:,cols.str.endswith('_traincounts')] = \
        total_adata.obsm[key_added].loc[:,cols.str.endswith('_traincounts')].fillna(0.0)
        total_adata.obsm[key_added][cols[cols.str.endswith('_traincounts')]] = \
        total_adata.obsm[key_added][cols[cols.str.endswith('_traincounts')]].astype(float)
    except Exception as e:
        logg.error(f'Error converting adata.obsm[{key_added}] contents to savable categories',exc_info=1)
    total_adata.obs_names = total_adata.obs['MMoCHi_obs_names']
    total_adata.obs.drop(columns='MMoCHi_obs_names', inplace=True)
    gc.collect()
    return total_adata,hierarchy

def _predict_probabilities_cutoff(X: sp.csr_matrix, clf,
                         probability_cutoff: float=0) -> Tuple[np.array, np.array, np.array]: 
    """
    Predict using the classifier, using the probability cutoff as described. Replace predictions with a ? if they do not match the cutoff.
    
    Parameters
    ----------
    X
        The matrix of cell features to predict on
    clf
        The classifier to use for predictions
    probability_cutoff: 
        The minimum confidence probability to tolerate before that cells' classification is removed. This feature is experimental.
        
    Returns
    -------
    full_probabilities: np.array
        n_classes x n_events array of the probabilities of each class for each event
    probability: np.array
        a 1 x n_events array of the probability used for each cell call. 
    predictions: np.array
        a 1 x n_events array of str, designating the class called for each group.
    """
    # predict_proba does not function with int64 encoded sparse matrix indices (caused by large datasets)
    if X.indices.dtype == np.int64:
        int_max = np.iinfo(np.int32).max
        batches = list()
        for i in range(0, X.shape[0], int_max // X.shape[1]): 
            batches.append(X[i:(i+int_max//X.shape[1])])
        for i in batches:
            i.indices.dtype = np.int32
    else:
        batches = [X]
        
    full_probabilities = list()
    for i in batches:
        probability = clf.predict_proba(i)
        full_probabilities.append(probability)

    full_probabilities = np.vstack(full_probabilities)
    
    probabilities = full_probabilities
    n_samples = probabilities.shape[0]
    class_type = clf.classes_.dtype    # all dtypes should be the same, so just take the first
    predictions = np.empty(n_samples, dtype=class_type)

    for k in range(n_samples):
        predictions[k] = clf.classes_[np.argmax(probabilities[k])]

    probabilities = np.max(probabilities,axis=1)

    for i,prob in enumerate(probabilities):
        if prob < probability_cutoff:
            predictions[i] = "?"
    return full_probabilities, probabilities, predictions

def hc_threshold(adata: anndata.AnnData, hierarchy: hierarchy.Hierarchy,
                 level: str, data_key: Optional[str]=utils.DATA_KEY, 
                 plot: bool=False, reference: Optional[str]=None,
                 batch: Optional[str]=None, log_reference: bool=True) -> pd.DataFrame:
    """
    Performs high-confidence thresholding using the subset definitions defined in one level of a MMoCHi hierarchy. 
    Note: This function relies on `mmc.utils.get_data()` for extracting expression information from the adata object. Thus, for each marker it will first search for matching proteins, then genes, and finally the .obs for matching feature names. If '_gex' or '_obs' are appended to the marker name, this priority order will be bypassed. 

    
    Parameters
    ----------
    adata:
        Object containing gene expression data, and expression data for modalities for every data key in .obsm
    hierarchy: Hierarchy object 
        Specifies one or more classification levels and subsets.
    level
        Name of a classification level in the hierarchy to threshold and high-confidence threshold on
    data_key
        Key in adata.obsm to be used for high-confidence thresholding
    plot
        Passed to thresholding.threshold, whether to display histograms of thresholds
    reference: 
        Column in the adata.obs to be used for logged debug info to reveal how the hc_threshold function performs       
    batch
        Name of a column in adata.obs that corresponds to a batch for use in the classifier   
    log_reference
        If true, logs high-confidence thresholding information
        
    Returns
    -------
    templin
        Dataframe of high-confidence calls for which classification an event falls into for the specified level of the hierarchy
    """
    # Get necessary information
    markers,values = hierarchy.classification_markers(level)
    
    # Setup templin
    templin = pd.DataFrame(index=adata.obs_names)
    temp = pd.DataFrame(index=adata.obs_names)
    for marker in markers:
        threshold = None
        try:
            interactive, preset_threshold = hierarchy.get_threshold_info(marker,level,batch)
        except:
            interactive, preset_threshold = True, None
        try:
            threshold, temp[marker] = thresholding.threshold(marker,adata,data_key,plot=plot,interactive=interactive,
                                                        preset_threshold=preset_threshold,run=True)
            
        except:
            # If thresholding doesn't work, attempt to find it anywhere else
            temp[marker] = utils.get_data(adata,marker,data_key)
            
        if not threshold is None:
            hierarchy.set_threshold(marker,threshold, False, level, batch)
    # Set up and excecute hc clauses
    clauses = []
    for label,truths in zip(values.keys(),values.values()):
        for single_truths in truths:
            clause = ""
            for mark, truth in zip(markers,single_truths):
                if truth != "any":
                    clause = clause + "(temp[\"" + mark + "\"] == \"" + truth + "\") & "
            clause = clause[:-3]
            if not reference is None:
                temp['mask'] = False
                exec("temp['mask'] = ("+clause+".values)",locals())
                if log_reference:
                    logg.debug((label,clause))
                    logg.debug('\n{}'.format(adata.obs[temp['mask'].values][reference].value_counts().to_string()))
                else:
                    logg.print([label,clause])
                    logg.print(adata.obs[temp['mask'].values][reference].value_counts())
                del temp['mask']
            clause = "temp.loc["+clause+",'names'] = \'"+label+"\'"
            clauses.append(clause)
    temp['names'] = "?"
    for clause in clauses:
        exec(clause,locals())
    
    templin[level+'_hc'] = temp['names']

    if not reference is None:
        df = pd.crosstab(templin[level+'_hc'],adata.obs[reference]).T
        if log_reference:
            logg.debug('high-confidence populations calculated...')
            logg.debug(f"\n{templin[level+'_hc'].value_counts()}")
            logg.debug('\n{}'.format(df.to_string()))
        else:
            logg.print('high-confidence populations calculated...')
            logg.print(templin[level+'_hc'].value_counts())
            logg.print(df)
    
    return templin
    
def _get_train_data(X: sp.csr_matrix, templin: pd.DataFrame, 
                    level: str, resample_method: Optional[str], 
                    max_training: Optional[int], in_danger_noise_checker: Union[bool,str], 
                    features_used: List[str]) -> Tuple[sp.csr_matrix,pd.DataFrame, np.array]:
    """
    Subsets sparse matrix of features into training and testing data.  If SMOTE or its variations are used, _traincounts will refer to rebalanced training counts and will not account for generated training events.
    
    Parameters
    ----------
    X: sparse matrix
        Subset of matrix of features to subset into training data
    templin: pandas DataFrame
        Dataframe of high-confidence calls for which classification an event falls into for the specified level of the hierarchy
    level: str
        Name of a classification level in the hierarchy to threshold and hc_threshold on
    resample_method: optional, str
        Method used to resample high-confidence data prior to training, passed to _balance_training_classes. "Oversample" or None. Other methods are experimental. See _balance_training_classes and the imblearn package for details. 
    max_training: optional, int
        If amount of training data is greater than max_training, a subset of the training data is randomly chosen to be used
    in_danger_noise_checker: bool or str
        Whether to check for in danger or noise values in the dataset
    features_used: listlike of str
        List of features that were checked/used in the reduction process 

    Returns
    ------- 
    X: sparse matrix
        Resampled feature matrix for use in training
    y: np.array()
        Resampled classifications for use in training, corresponding to rows in X
    traincounts_final:
        array of length templin with the number of times each event is used in training
    """
    logg.print('Choosing training data...') # Randomly select high-confidence dataset for training purposes...80% train, 20% held out
    # creating an array to later fill with traincounts by applying masks in reverse order of y subsetting
    traincounts_final = np.zeros(len(templin))
    templin[level+'_hold_out'] = True
    mask = np.random.choice([False,False,False,False,True],size=sum(templin[level+'_hc'] != '?'))
    templin.loc[templin[level+'_hc'] != '?',level+'_hold_out'] = mask
    logg.debug(f"Number of events {len(templin)}")
    logg.debug(f"X.shape: {X.shape}")
    X = X[~templin[level+"_hold_out"].to_numpy()]
    y = templin[~templin[level+'_hold_out']][level+'_hc'].values
    traincounts_mid = np.zeros(len(y))
    idx_remove_list=[]
    n_idxs = np.array(range(0,len(y)))
    max_training_per_subset = int(max_training / len(set(y)))
    for i in set(y):
        # Remove from classes above max_training size
        total = sum(y==i)
        if total > (max_training_per_subset):
            idxs = n_idxs[y==i]
            idx_remove = np.random.choice(idxs,total-max_training_per_subset,replace=False)
            idx_remove_list.append(idx_remove)
    try:
        idx_remove = np.concatenate(idx_remove_list)
        idxs_remove_mask = np.array([j in idx_remove for j in n_idxs])
    except:
        idxs_remove_mask = np.array([False] * len(y))
        logg.debug(f'None removed at {level}')
    y = y[~idxs_remove_mask]
    X = X[~idxs_remove_mask]
    logg.info(f"{len(y)} real cells in training set...")
    logg.print('Resampling...')
    X,y, tcounts = _balance_training_classes(X,y,features_used, resample_method,max_training=max_training,in_danger_noise_checker=in_danger_noise_checker) 
    traincounts_mid[~idxs_remove_mask] = tcounts
    traincounts_final[~templin[level+'_hold_out']] = traincounts_mid 
    logg.debug(f'Selected for training: {len(y)}, with {np.unique(y, return_counts=True)}')
    logg.print(f"Training with {len(y)} events after {resample_method} resampling...")
    assert len(set(y)) > 1, 'Must have more than 1 class to classify...'

    return X, y, traincounts_final

def _in_danger_noise(nn_estimator, samples: List[Tuple[anndata.AnnData, sp.csr_matrix]],
                     y: Union[anndata.AnnData, sp.csr_matrix]) -> np.array:
    """
    Estimate if events in a sample are in danger or noise. This calculation is performed by clustering and performing nearest neighbors on samples to identify events with no niehgbors in agreement (noise) or partial neighbor agreement (in danger).
    Adapted from BorderlineSMOTE and SVMSMOTE of the imblearn package.
    Note that this function relies on scanpy's implementations of HVGs and PCA, but there may be a better approach computationally for this.
    
    Parameters
    ----------
    nn_estimator : estimator object
        An estimator that inherits from
        :class:`~sklearn.neighbors.base.KNeighborsMixin` use to determine
        if a sample is in danger/noise.
    samples : {array-like, sparse matrix} of shape (n_samples, n_features)
        The samples to check if either they are in danger or not
    target_class : int or str
        The target corresponding class being over-sampled
    y : array-like of shape (n_samples,)
        The true label in order to check the neighbor labels
        
    Returns
    -------
    items: np.Array
        for each data point in the array has a value of
            0 if neither in danger nor in noise
            1 if just in danger
            2 if just noise
            3 if both
    """
    logg.debug('kneighbors for in danger noise')
    x = nn_estimator.kneighbors(samples, return_distance=False)[:, 1:]
    items = np.zeros(samples.shape[0])
    for target_class in set(y):
        nn_label = (y[x] != target_class).astype(int)
        n_maj = np.sum(nn_label, axis=1) #number of different class adjacent

        # Samples are in danger for m/2 <= m' < m
        danger = (n_maj >= (nn_estimator.n_neighbors - 1) / 2) & (n_maj <= nn_estimator.n_neighbors - 4) & (y == target_class)
        # Samples are noise for m = m'
        noise = ((n_maj >= nn_estimator.n_neighbors - 2)  & (y == target_class)).astype(int) * 2
        items = items+danger+noise

    # Targets can also be labeled "in danger" if they are clustered in a small cluster that's not well-represented
    leiden = nn_estimator.cluster()
    crosstab = pd.crosstab(leiden,y)
    logg.debug(crosstab)
    crosstab = ((crosstab.div(crosstab.sum(axis=0),axis=1)) <= .05) & ((crosstab.div(crosstab.sum(axis=1),axis=0)) >= .75) & (crosstab >= 5)
    for i in crosstab:
        for j,k in zip(crosstab[i],crosstab[i].index):
            if j:
                items[(leiden == k) & (y == i) & (items==0)] = 2
    return items
    
class _HVF_PCA_Neighbors():
    def __init__(self,features_used,**kwargs):
        self.features_used = np.array(features_used)
        self.feature_modalities = np.array([x.split('_mod_')[-1] for x in features_used])
        self.modalities = Counter(self.feature_modalities)
        self.neighbors = sklearn.neighbors.NearestNeighbors(**kwargs)
        self.n_neighbors = self.neighbors.n_neighbors
        self.leiden = None
    def fit(self, X: Union[np.array, sp.csr_matrix]):
        """
        Runs PCA on the highly variable features of self
        
        Parameters
        ----------
        X: numpy array or sparse matrix 
            Self.X
        Returns
        -------
        self: _HVF_PCA_Neighbors object
            Modified with X_pca column in .obsm with PCA of the highly variable features
        """
        self.hvf = []
        for i in self.modalities:
            logg.debug(f'Computing HVFs for {i}')
            if self.modalities[i] < 5000:
                self.hvf.extend(self.features_used[self.feature_modalities == i])
            else:
                adata = anndata.AnnData(X=X[:,self.feature_modalities == i],dtype=np.float32, 
                                        var = pd.DataFrame(index=self.features_used[self.feature_modalities == i]))
                sc.pp.highly_variable_genes(adata,n_top_genes = 5000)
                self.hvf.extend(adata.var_names[adata.var.highly_variable].tolist())
        logg.debug('PCA for fit for _HVF_PCA_Neighbors')   
        adata = anndata.AnnData(X=X[:,[i in self.hvf for i in self.features_used]],dtype=np.float32)
        sc.pp.scale(adata,zero_center=False)
        sc.tl.pca(adata,n_comps=15)
        self.pca = adata.obsm['X_pca']
        X = self.pca[:,0:15].copy()
        del adata
        self.neighbors = self.neighbors.fit(X)
        self.leiden = None
        return self
    def kneighbors(self,*args,**kwargs):
        return self.neighbors.kneighbors(**kwargs)
    def kneighbors_graph(self,*args,**kwargs):
        return self.neighbors.kneighbors_graph(**kwargs)
    def cluster(self):
        """
        TODO implement neighbors and leiden calculations that do not depend on scanpy"""
        if self.leiden is None:
            adata = anndata.AnnData(X=self.pca,dtype=np.float32)
            adata.obsm['X_pca'] = self.pca
            logg.debug('neighbors for in_danger_noise clustering')
            sc.pp.neighbors(adata,n_pcs=15)
            logg.debug('clustering for in_danger_noise')
            # skip aall this and do adjacency=kneighbors_graph
            sc.tl.leiden(adata,resolution = .5,n_iterations=1)
            if len(adata.obs.leiden.unique()) > min(10,len(adata)/100):
                sc.tl.leiden(adata,resolution = .25,n_iterations=1)

            self.leiden = adata.obs.leiden.values
        return self.leiden
    def get_params(self,**kwargs):
        return self.neighbors.get_params(**kwargs)

def _borderline_balance(nn_estimator, X: sp.csr_matrix,
                       y: np.array, indices: np.array, min_cluster_size: float=0.01) -> Tuple[sp.csr_matrix, np.array, np.array]:
    """
    Uses nearest neighbors clustering to cluster data in y, and for clusters larger than the minimum cluster size, duplicates events in training data to equalize its representation.
    
    Parameters
    ----------
    nn_estimator
        nearest neighbor clustering estimator, must have cluster method
    X
        .X column that references the same AnnData object as Y
    y
        Data to cluster and balance
    indices
        indices of initial X and y passed to _balance_training_classes for keeping track of training counts
    min_cluster_size
        Minimum proportion of all classes that a cluster must contain
    
    Returns
    -------
    X: scipy.sparse.csr_matrix()
        .X column of clustered data
    y: np.Array
        Clustered data
    indices
        subsetted indices
    """
    leiden = nn_estimator.cluster()
    crosstab = pd.crosstab(leiden,y)
    mins = crosstab.div(crosstab.sum(axis=0)).min(axis=1)
    mults = crosstab.div(crosstab.sum(axis=0)).div(mins,axis=0)
    leidens = []
    for i in set(leiden):
        # If a cluster makes up at least 1% of all classes...
        if mins[i] > min_cluster_size:
            # Equalize the representation of that cluster across the board
            for j in mults.columns: 
                mult = int(np.ceil(mults.at[i,j]-1))
                if mult > 0:
                    mask = (leiden==i)&(y==j)
                    X = sp.vstack([X]+[X[mask]]*mult,format='csr')        
                    y = np.hstack([y]+[y[mask]]*mult)
                    indices = np.hstack([indices]+[indices[mask]]*mult)
                    leiden = np.hstack([leiden]+[leiden[mask]]*mult)
            leidens.append(i)
            
    mask = [i in leidens for i in leiden]
    if (sum(mask) < (.5*len(mask))) & (sum(mask)>0):
        mult = int(np.ceil((.5*len(mask)/sum(mask))))
        X = sp.vstack([X]+[X[mask]]*mult,format='csr')        
        y = np.hstack([y]+[y[mask]]*mult)
        indices = np.hstack([indices]+[indices[mask]]*mult)
        leiden = np.hstack([leiden]+[leiden[mask]]*mult)
    return X, y, indices
            
def _balance_training_classes(X: sp.csr_matrix, y: np.array,
                              features_used: List[str], method: Optional[str]=None,
                              max_training: Optional[int]=None, 
                              in_danger_noise_checker: Union[bool,str]=True) -> Tuple[sp.csr_matrix, np.array, np.array]:
    """
    Methods to balance classes, using the imblearn package. Consult their documentation for more details
    
    Parameters
    ----------
    X: sparse matrix
         .X column that references the same AnnData object as y
    y: np.Array
        Data to sample and balance
    features_used: Listlike of str
        List of features that were checked/used in the reduction process 
    method: str, possible methods: "undersample", 'resample', 'oversample', 'SMOTE', 'KMeansSMOTE', 'BorderlineSMOTE' ,None
        Method used to balance training classes (over or undersampling). See imblearn for details. If SMOTE or its variations are used, _traincounts will refer to rebalanced training counts and will not account for generated training events.
    max_training: int
        Maximum length of y data to be used in training, if len(y) > max_training will randomly subsample x and y
    in_danger_noise_checker: bool or str, if str can include "danger", "noise", or both
        Whether to remove noise data and boost in danger data points
        True is equivalent to "danger and noise"
        
    Returns
    -------
    X: scipy.sparse.csr_matrix
        Resampled feature matrix, with balanced classes
    y: np.array()
        Resampled classifications for use in training, corresponding to rows in X
    traincounts: 
        len(x) with the number of times each value is sampled
    """
    full_indices = np.array(range(len(y)))
    full_len = len(y)
    if in_danger_noise_checker or method in ['SMOTE','BorderlineSMOTE','KMeansSMOTE']:
        logg.debug('begining HVF PCA NEIGHBORS')
        k_neighbors = _HVF_PCA_Neighbors(features_used,n_jobs=-1,n_neighbors=10).fit(X)

    if in_danger_noise_checker:
        if type(in_danger_noise_checker) is bool:
            in_danger_noise_checker = 'danger and noise'
        logg.debug('begining in_danger_noise')
        danger_or_noise = _in_danger_noise(k_neighbors,X,y)
        logg.info(f'Found: {sum(danger_or_noise==2)} noise and {sum(danger_or_noise==1)} in danger of {len(danger_or_noise)} events.')
        logg.debug(pd.crosstab(danger_or_noise,y))
        if 'noise' in in_danger_noise_checker:
            if set(y[danger_or_noise != 2]) == set(y):
                full_indices = full_indices[danger_or_noise != 2]
                y = y[danger_or_noise != 2]
                X = X[danger_or_noise != 2]
                if not k_neighbors.leiden is None:
                    k_neighbors.leiden = k_neighbors.leiden[danger_or_noise != 2]
                danger_or_noise = danger_or_noise[danger_or_noise != 2]
            else:
                logg.info('Could not remove noise, as it would remove a category')
        else:
            logg.print('Not looking for noise')
        if 'danger' in in_danger_noise_checker:
            if sum(danger_or_noise) > 0:
                X = sp.vstack([X]+[X[danger_or_noise==1]]*4,format='csr')     
                y = np.hstack([y]+[y[danger_or_noise==1]]*4)
                full_indices = np.hstack([full_indices]+[full_indices[danger_or_noise==1]]*4)
                if not k_neighbors.leiden is None:
                    k_neighbors.leiden = np.hstack([k_neighbors.leiden]+[k_neighbors.leiden[danger_or_noise==1]]*4)
            else:
                logg.debug('No "in danger" found')
        else:
            logg.debug('Not looking for danger')
            
        if 'rebalance' in in_danger_noise_checker:
            X, y, full_indices = _borderline_balance(k_neighbors,X,y, full_indices)
    
    
    n_tcounts = len(y)
    if method == "oversample":
        sample = imblearn.over_sampling.RandomOverSampler()
    elif method == "resample":
        sample = imblearn.under_sampling.OneSidedSelection()
    elif method == "undersample":
        sample = imblearn.under_sampling.RandomUnderSampler()
    elif method == "SMOTE":
        sample = imblearn.over_sampling.SMOTE(k_neighbors=k_neighbors,n_jobs=-1)
    elif method == "BorderlineSMOTE":
        sample = imblearn.over_sampling.BorderlineSMOTE(k_neighbors=k_neighbors,n_jobs=-1)
    elif method == "KMeansSMOTE":
        sample = imblearn.over_sampling.KMeansSMOTE(k_neighbors=k_neighbors,n_jobs=-1) 
    else:
        logg.info('No sampling method chosen')
    if method in ["undersample",'resample','oversample','SMOTE','KMeansSMOTE','BorderlineSMOTE']:
        X, y = sample.fit_resample(X, y)
        if 'SMOTE' in method:
            training_indices = np.array(range(n_tcounts))
        else:
            training_indices = sample.sample_indices_
    else:
        training_indices = np.array(range(n_tcounts))
    
    
    if not max_training is None and len(y) > max_training: 
        X,y, training_indices = sklearn.utils.resample(X,y, training_indices, replace=False,n_samples=max_training)

    training_indices = full_indices[training_indices]
    tcounts = np.bincount(training_indices.astype(int),minlength=full_len)
    tcounts = tcounts
    return X, y, tcounts

def terminal_names(adata: anndata.AnnData, obs_column: str='classification',
                   confidence_column: str='certainty', conf_drops: Union[float, List[float]]=None,
                   key_added: str='lin', voting_reference: str = None, 
                   majority_min: float=.9) -> anndata.AnnData:
    """
    Create a column in the .obs featuring the most specific classification/subset for each event. Uses horizonal concatenation (rightward has priority) of classifications in the adata.obsm[key_added]
    
    Parameters
    ----------
    adata: AnnData object 
        Object containing gene expression data, and expression data for modalities for every data key in .obsm
    obs_column: str
        Column in the .obs to add to the AnnData containing the classification 
    confidence_column: str
        Column in the .obs to add to the AnnData containing the percent of trees in agreement with the classification 
    conf_drops: float or listlike of float
        Maximum percent of confidence-loss between levels of the hierarchy before an event should be terminated at an earlier subset
    key_added: str
        Key in the adata object that contains the classification results
    voting_reference: str
        Column in the .obs corresponding to leiden clusters or overclusters for use in the voting function (if majority voting)
    majority_min: float [0,1]
        Proportion of events within an individual voting group that must agree for that voted identity to override all members of that voting group
    
    Returns 
    -------
    AnnData: AnnData object
        Contains a new column for classifications and potentially a column for confidences
    """
    # ID classification columns
    class_names = [x for x in adata.obsm[key_added].columns if "_class" in x][1:]
    names_matrix = adata.obsm[key_added].loc[:,class_names].astype(str)
    names_matrix.columns = [i.split('_class')[0] for i in names_matrix.columns]
    # ID confidence columns
    conf_names = [x for x in adata.obsm[key_added].columns if "_probability" in x]
    # Match column names
    conf_matrix = pd.DataFrame(columns = names_matrix.columns, index = names_matrix.index)
    probability_matrix = adata.obsm[key_added].loc[:,conf_names]
    probability_matrix.columns = [i.split('_probability')[0] for i in probability_matrix.columns]
    conf_matrix = conf_matrix.fillna(probability_matrix)
    conf_matrix = conf_matrix.fillna(1)
    conf_matrix[(names_matrix.isna()) | (names_matrix == 'nan')] = np.nan
    
    if not confidence_column is None:
        conf_matrix.columns = names_matrix.columns
    
    if conf_drops is None:
        conf_drops = -1

    if not isinstance(conf_drops, list):
        conf_drops = [conf_drops]*len(names_matrix.columns)
    
    labels = names_matrix.iloc[:,0].copy()
    labels.name = None
    confidences = conf_matrix.iloc[:,0].copy()
    confidences.name = None
    blacklisted = conf_matrix.iloc[:,0].copy() < 0
    blacklisted.name = None

    for col, conf_drop in zip(names_matrix.columns,conf_drops):
        mask = conf_matrix.loc[:,col].copy()
        mask.loc[blacklisted] = -1

        mask.name = None
        labels.loc[mask>conf_drop] = names_matrix.loc[mask>conf_drop,col]
        confidences.loc[mask>conf_drop] = conf_matrix.loc[mask>conf_drop,col]
        blacklisted.loc[mask<conf_drop] = True

    adata.obs[obs_column] = labels
    if not confidence_column is None:
        adata.obs[confidence_column] = confidences
    adata.obs[obs_column] = adata.obs[obs_column].astype('category')
    if voting_reference is None:
        return adata
    else:
        def percent_mode(x):
            mode = x.value_counts().index[0]
            return sum(x == mode)/len(x)
        simple_majority = adata.obs.groupby(voting_reference).agg({obs_column:lambda x:x.value_counts().index[0]})[obs_column]
        majority = simple_majority[(adata.obs.groupby(voting_reference).agg({obs_column:percent_mode}) > majority_min)[obs_column]]
        adata.obs['simple_majority_voting'] = adata.obs[voting_reference].map(simple_majority)
        adata.obs['majority_voting'] = adata.obs[voting_reference].map(majority).fillna(adata.obs[obs_column])
        return adata
    
def identify_group_markers(adata: anndata.AnnData, group1: Union[str, Iterable[str]],
                           group2: Union[str, Iterable[str]]=[],
                           batch_val: Optional[Union[str, Iterable[str]]]=None, 
                           batch: str=utils.BATCH_KEY,
                           reference: str='leiden', key_added: str='groups',
                           filtered: bool=False, plot: bool=True, 
                           min_fold_change: int=2, min_in_group_fraction: float=.5,
                           max_out_group_fraction: float=0.25, n_ups: int=29,
                           n_downs: int=30, use_raw: bool=False, return_df: bool=True) -> pd.DataFrame:
    """
    Calculates differentially expressed genes between two provided groups (specified in the .obs).
    Built on scanpy.tl.rank_genes_groups
    Use https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html for further details and documentation.
    
    Parameters
    ----------
    adata: AnnData object 
        Log normalized data to run differential expresssion on. Contains group 1 and 2 in .obs. Must have X_umap to plot.
    group1: Union[str, Iterable[str]]
        Subset of groups to compare on, see sc.tl.rank_genes_groups for more detailed documentation
    group2: Union[str, Iterable[str]]
        Second subset of groups to compare on, see sc.tl.rank_genes_groups for more detailed documentation
    batch_val: str or list of str
        Batch or batches used to limit the events that differential expression is run on (in adata.obs[batch])
    batch: str
        Name of a column in adata.obs that corresponds to a batch for use in the classifier   
    reference: str
        .obs tag to find and split groups 1 and 2 on
    key_added: str
        Key in the adata object that contains the classification results
    filtered: bool
        Uses sc.tl.filter_rank_genes_groups to filter out genes base on log fold changes and genes inside and outside groups
    plot: bool
        Whether to plot a UMAP of the events colored by their differential expression, using X_umap. 
    min_fold_change: int
        If filtered, passed to sc.tl.filter_rank_genes_groups
    min_in_group_fraction: float [0,1]
        If filtered, passed to sc.tl.filter_rank_genes_groups
    max_out_group_fraction: float [0,1]
        If filtered, passed to sc.tl.filter_rank_genes_groups
    n_ups: int
        If plot, includes best n_ups values
    n_downs: int
        If plot, includes worst n_downs values
    use_raw: bool
        Whether to pull expression for DE analysis from the .raw
    return_df: bool
        Whether to return a df of genes ranked in their ability to characterize group 1
    Returns
    -------
    df: pandas df 
        If return_df, returns pandas df of genes ranked in their ability to characterize group 1 
    """
    adata.obs[key_added] = 'outgroup'
    group1 = group1 if isinstance(group1,list) else [group1]
    group2 = group2 if isinstance(group2,list) else [group2]
    if not batch_val is None:
        batch_val = batch_val if isinstance(batch_val,list) else [batch_val]

    if reference == 'leiden':
        group1, group2 = [str(i) for i in group1], [str(i) for i in group2]
        
    group1_mask = adata.obs[reference].isin(group1)
    
    if group2 == []:
        group2_mask = ~adata.obs[reference].isin(group1)
    else:
        group2_mask = adata.obs[reference].isin(group2)
    
    if (batch_val is None) or (batch is None):
        batch_mask = [True]*len(adata)    
    else:
        batch_mask = adata.obs[batch].isin(batch_val)
        group1_mask, group2_mask = (group1_mask & batch_mask), (group2_mask & batch_mask)
        
        
    adata.obs.loc[group1_mask,key_added] = 'group1'
    adata.obs.loc[group2_mask,key_added] = 'group2'
    adata.obs[key_added] = adata.obs[key_added].astype('category')
    
    sc.tl.rank_genes_groups(adata, groupby=key_added, groups=['group1'], reference='group2', use_raw=use_raw)
    
    if filtered:
        try:
            sc.tl.filter_rank_genes_groups(adata, key='rank_genes_groups', key_added='rank_genes_groups', groupby='groups', compare_abs=True, 
                                       min_fold_change=min_fold_change, min_in_group_fraction=min_in_group_fraction, 
                                       max_out_group_fraction=max_out_group_fraction)
        except:
            sc.tl.filter_rank_genes_groups(adata, key='rank_genes_groups', key_added='rank_genes_groups', groupby='groups', 
                                       min_fold_change=min_fold_change, min_in_group_fraction=min_in_group_fraction, 
                                       max_out_group_fraction=max_out_group_fraction)
        
    df=sc.get.rank_genes_groups_df(adata, group='group1')
    df= df[~df.names.isna()]
    if plot:
        if len(df) > n_ups+n_downs:
            ups,downs = df.head(n_ups).names.to_list(), df.tail(n_downs).names.to_list()
        else:
            ups,downs = df.names.to_list(), []
        sc.pl.umap(adata[batch_mask].copy(), color=[key_added]+ups+downs, cmap='jet', s=10,ncols=6)
        gc.collect()
    
    if return_df:
        return df
    return