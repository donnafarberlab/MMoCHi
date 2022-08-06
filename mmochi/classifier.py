import numpy as np
import pandas as pd
import sklearn 
import sklearn.ensemble
import sklearn.neighbors

import imblearn
import random
import gc
from collections import Counter
import scipy.sparse as sp
import scipy.stats as stats
import statistics 
import scanpy as sc
import anndata

from . import utils
from .logger import logg
from . import thresholding 

def classifier_setup(adata, data_key=None, reduce_features_min_cells=0, features_limit=None):
    '''Setup before training a classifier. This can be run before the classifier (to reduce runtime of the classifier function in a parameter optimization loop),
    or is automatically run when training a classifier. It concatenates the .X and any data_key in the .obsm, then performs feature reduction (if reduce_features_min_cells > 0).
    Next, features can be limited by an external feature set.
    Then, it sorts the resulting feature_names (the columns from the .X and .obsm[data_key]) and csr.matrix, alphabetically, to make the feature order reproducible across runs. 
    If defined, feature limits can be performed so that you can match the expected features of the hierarchy
    TODO, add in NA checks
    
    adata: anndata object
    data_key: str, key in adata.obsm
    reduce_features_min_cells: int, passed to _reduce_features
    features_limit: listlike of str
    '''   
    print('Setting up...')
    
    if data_key is None or features_limit == 'GEX':
        print('Using .X only')
        X = sp.csr_matrix(adata.X)
        features_used = list(adata.var_names)
    else:
        print('Using .X and',data_key)
        X = sp.hstack([adata.X,adata.obsm[data_key]],format='csr')
        features_used = list(adata.var_names) + list(adata.obsm[data_key].columns)
        
    print("Reducing features...")
    
    X, features_used = _reduce_features(X,features_used, min_cells=reduce_features_min_cells)
    
    if not features_limit is None and not features_limit == 'GEX':
        logg.warn("Limiting features...")
        feature_mask = np.isin(features_used, feature_names)
        logg.debug(sum(feature_mask))
        features_used = np.array(features_used)[feature_mask]
        logg.debug(X.shape)
        X=X[:,np.flatnonzero(feature_mask)]
        logg.debug(X.shape)
    features_used = np.array(features_used)
    sort = np.argsort(features_used)
    X = X[:,sort]
    logg.debug(features_used)
    features_used = features_used[sort]
    
    return X,list(features_used)

def _reduce_features(X,feature_names, min_cells=50):
    '''Reduce features that only vary in fewer than min_cells_threshold cells.
    This is a very simplistic feature reduction method. Planned additions include batch based methods and other options
    
    X: 2D sp.csr_matrix 
    feature_names: listlike of feature names 
    min_cells: minimum cells threshold. Any feature expressed in fewer cells will be excluded.
    '''
    if min_cells > 0:
        num_cells = X.getnnz(axis=0)
        X=X[:,num_cells>min_cells]
        feature_names = np.array(feature_names)[num_cells>min_cells]
    else:
        print("No min_cells feature reduction...")
    return X, feature_names

def classify(adata,hierarchy,key_added='lin',data_key='protein',x_data_key=None,
             batch_key=None,retrain=False,plot_thresholds=False,
             reduce_features_min_cells=25, allow_spike_ins=True,
             score=False, proba_suffix='_proba', # for UNS. I should honestly fix all of these...
             resample_method='oversample',X=None,features_used=None,in_danger_noise_checker=True,
             probability_cutoff=0,max_training=20000,features_limit = None, skip_to = None, end_at = None,
             n_estimators=100,reference=None, clf_kwargs={}):
    ''' Build a classification from a hierarchy. If the hierarchy has not been trained yet, this function trains the hierarchy AND 
    predicts the dataset. If the hierarchy has been trained, this classifier will skip retraining unless explicitly asked for. 
    adata: anndata object
    hierarchy: Hierarchy object specifying one or more classification levels and subsets.
    key_added: str, key in the .obsm to store a dataframe of the classifier's results
    proba_suffix: str, If defined, probability outputs of predict_proba will be saved to the .uns[level+proba_suffix]
    
    retrain: bool, If the classification level of a hierarchy is already defined, whether to replace that classifier with a retrained one.
               If False, no new classifier would be trained, and that classifier is used for prediction of the data.
    
    reference: str, column in the adata.obs to be used for logged debug info to reveal how the ground truth function performs
    
    data_key: str, key in adata.obsm to be used for ground truth thresholding
    x_data_key: str, key in adata.obsm, used for the creation of the training data, if undefined, defaults to data_key
    batch_key: str, name of a column in adata.obs that corresponds to a batch for use in the classifier
                   
    plot_thresholds: bool, Whether to display thresholds when calculating ground truth
    score: bool, Whether to plot confusion matricies after each training and classification level
    
    reduce_features_min_cells: int, remove features that are expressed in fewer than this number, passed to _reduce_features
                                    feature reduction can be a very powerful tool for improving classifier performance
    features_limit: str, listlike of str, either "GEX" or a list of str specifying features to limit to
    
    X, features_used: sp.csr_matrix and listlike of str, Optional setup data to circumvent the initial classifier setup functions
                                                         The presence of both of these overrides any predefined feature reduction techniques.
    
    probability_cutoff: float, Must be between 0 and 1. The minimum fraction of trees that must agree on a call for that event to be labelled 
                               by the classififer.     

    resample_method: str, method to resample ground truth data prior to training, passed to _balance_training_classes
    max_training: int, maxiumum number of events to use in training. Anything more, and training data will be subsampled to meet this number
                       a high max_training increases the computational time it takes to train the model, but may improve reliability
    n_estimators
    clf_kwargs: dict, default kwargs (can be overridden by subset kwargs on the hierarchy)
    
    allow_spike_ins: bool, Whether to allow spike ins when doing batched data. Spike ins are performed if a subset is below the minimum threshold
                           within an individual batch, but not the overall dataset. If False, errors may occur if there are no cells of a particular subset
                           within a batch. Warning, spike ins currently do not respect held out data, meaning there will be much fewer than 20% held out data if 
                           spike ins are frequent. # TODO, FIX
    
    skip_to, end_at: str, names of hierarchy classification levels to start at or prematurely end at. This is extremely useful when
                          debugging particular levels of the classifier. The order of levels corresponds to the order they are defined
                          in the hierarchy. Requires that the anndata object has a predefined .obsm[key_added]. If skip_to is defined, 
                          it replaces only the columns that occur at and beyond that level.
        
    returns: adata with a .obsm[key_added] containing columns corresponding to the classification ground truth ("_gt"), the training data ("_train"), 
             the prediction ("_class") for each level. If proba_suffix is defined, also contains "_proba" keys in the .uns corresponding to the prediction 
             probabilities of each class for each event.
             Also returns a trained hierarchy with classifiers built into it and a record of thresholds and other settings used.
             
             Recommended to follow the use of this function with cl.terminal_names().
    '''
    if x_data_key is None:
        x_data_key = data_key
    setup_complete = (not X is None) and (not features_used is None)
    if not setup_complete:
        X,features_used = classifier_setup(adata,x_data_key,reduce_features_min_cells,features_limit=features_limit)
    print(f"Set up complete.\nUsing {len(features_used)} features")
    
    classifications = hierarchy.get_classifications()
    if skip_to is None:
        adata.obsm[key_added] = pd.DataFrame(index=adata.obs.index)
        adata.obsm[key_added]['All_class'] = 'All'
    else:
        assert skip_to in classifications, 'Attempting to skip to an invalid classificaiton level'
        classifications = classifications[classifications.index(skip_to):]
        parent,gparent = hierarchy.classification_parents(skip_to)
        assert gparent+"_class" in adata.obsm[key_added].columns, 'Attempting to skip to a level beneath an unclassified level'
        logg.warning(f'Skipped to classification on {skip_to}')
        for classification in classifications:
            adata.obsm[key_added].drop([classification+"_class",classification+"_gt",
                                        classification+"_train",classification+"_proba",classification+"_trains",classification+"_tcounts"],axis=1, inplace=True, errors='ignore')
    if not end_at is None:
        classifications = classifications[:classifications.index(end_at)+1]
            
    total_adata, total_X = adata, X
    total_adata_obs_names = total_adata.obs_names
    total_adata.obs_names = range(0,len(total_adata))    
    total_adata.obs_names = total_adata.obs_names.astype(str)

    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    
    for level in classifications:
        try:
            parent,gparent = hierarchy.classification_parents(level)
            assert gparent+"_class" in total_adata.obsm[key_added].columns, f"Cannot subset on {parent} to train {level} as {gparent} has not been classified"
            subset_mask = total_adata.obsm[key_added][gparent+'_class'] == parent
            assert sum(subset_mask) > 0, f"No events were classified as {parent} in {gparent}_class"
            subset_adata, subset_X = total_adata[subset_mask],total_X[subset_mask.values] 
            print("Data subsetted on", parent,'in', gparent) 
            
            if retrain or not hierarchy.has_clf(level):
                print(f'Running ground truths for {level}...')         
                gt_dfs, gt_batches = [], []
                for batch_mask, batch in zip(batch_masks, batches):
                    adata = subset_adata[batch_mask[subset_mask]]
                    if not batch_key is None:
                        logg.info(f'Running ground truths in {batch}')
                    if len(adata) > 0:
                        df = ground_truth(adata,hierarchy,level,data_key,reference=reference,
                                          batch=batch,plot=plot_thresholds)
                        gt_dfs.append(df)
                        gt_batches.append(batch)
                    else:
                        logg.info(f'No events in subsetted adata, skipping ground truth')
                    del adata
                    gc.collect()

                if hierarchy.is_cutoff(level):
                    print(f'Performing cutoff for {level}...')
                    df = pd.concat(gt_dfs)
                    df[level+'_class'] = df[level+'_gt']
                    df.loc[df[level+'_class'] == '?',level+'_class'] = np.nan
                else:
                    print(f'Preparing training data for {level}...')
                    gt_df = pd.concat(gt_dfs)
                    min_events = hierarchy.get_min_events(level)
                    
                    if min_events != 0:
                        logg.info(f'Checking subsets for minimum events...')
                        if min_events < 1:
                            min_events = np.ceil(min_events * len(subset_adata))
                            
                        subsets = hierarchy.subsets_info(level).keys()
                        vals = gt_df[level+'_gt'].value_counts()
                        vals.drop('?',errors='ignore',inplace=True)
                        logg.debug(f'Values before min_events ({min_events}) thresholding\n{vals.to_string()}')
                        logg.debug(f'len(vals) = {len(vals)}, len(subsets) = {len(subsets)}')
                        if min(vals) < min_events or len(vals) < len(subsets):
                            logg.debug(f'Must remove at least one subset')
                            removed_subsets = list(set(subsets)-set(vals.index[vals>=min_events]))
                            logg.warning(f'Removing {removed_subsets} cells as they are below {min_events} events')
                            gt_df.loc[gt_df[level+'_gt'].isin(removed_subsets),level+'_gt'] = '?'
                            subsets = list(set(subsets)-set(removed_subsets))
                        assert len(gt_df[gt_df[level+'_gt']!=0][level+'_gt'].unique()) > 1, 'Cannot train classifier with only one class'
                        
                        if not batch_key is None:
                            gt_barcodes_batches = []
                            for df, batch in zip(gt_dfs, gt_batches):
                                gt_barcodes_batch = list(df.index)
                                for subset in subsets:
                                    vals = df[level+'_gt'].value_counts()
                                    if not subset in vals.index or (vals.loc[subset] < min_events):
                                        if not subset in vals.index:
                                            events_needed = min_events
                                        else:
                                            events_needed = min_events - vals.loc[subset]
                                        if allow_spike_ins:
                                            logg.info(f'Spiking in {events_needed} of {subset} in {batch} to reach {min_events} events')
                                            gt_barcodes_batch += random.sample(list(gt_df[(gt_df[level+'_gt']==subset) & ~gt_df.index.isin(gt_barcodes_batch)].index), events_needed)
                                        else:
                                            assert False, f"{batch_key} of {batch} is missing subset {subset}. Cannot train without spike ins."
                                gt_barcodes_batches.append(gt_barcodes_batch)
                        else:
                            gt_barcodes_batches = [list(gt_df.index)]
                            
                    class_weight = hierarchy.get_class_weight(level)
                        
                    if class_weight == 'balanced' and not batch_key is None:
                        y_org = gt_df.loc[gt_df[level+'_gt']!='?',level+'_gt'].values
                        classes = np.unique(y_org)
                        class_weight = sklearn.utils.class_weight.compute_class_weight('balanced', classes=np.unique(y_org), y=y_org)
                        class_weight = dict(zip(classes,class_weight))
                        class_weight = {}
                        for subset in classes:                            
                            # class_weight[subset]= math.log((1/(sum(y_org==subset)/len(y_org))+1),100)
                            class_weight[subset] = (1/(sum(y_org==subset)/len(y_org))) ** (1./2)
                        logg.info(f'Manually made a balanced class_weight: {class_weight}')

                    clf_kwargss = _default_kwargs(clf_kwargs, {'n_jobs':-1, 'bootstrap':True, 'verbose':True})
                    clf_kwargss = _default_kwargs(hierarchy.get_kwargs(level).copy(),clf_kwargss)
                    logg.info(clf_kwargss)
                    print(f'Initializing classifier for {level}...')
                    clf = sklearn.ensemble.RandomForestClassifier(class_weight=class_weight,
                                                                  warm_start=(not batch_key is None),
                                                                  n_estimators=n_estimators,**clf_kwargs)
                    

                    dfs = []
                    train_dfs = []
                    for gt_barcodes_batch, batch in zip(gt_barcodes_batches, gt_batches):
                        logg.info(f'Training using {batch}...')
                        df = gt_df.loc[gt_barcodes_batch].copy()
                        mapper = dict(zip(subset_adata.obs_names,range(0,len(subset_adata.obs_names))))
                        gt_barcodes_batch = [mapper[i] for i in gt_barcodes_batch]
                        X,y = _get_train_data(subset_X[gt_barcodes_batch],df,level,resample_method,max_training,in_danger_noise_checker)
                        train_df = _map_training(subset_X[gt_barcodes_batch], X, df.index)
                        train_dfs.append(train_df)
                        clf.fit(X, y)
                        dfs.append(df)
                        clf.n_estimators += n_estimators
                        del X, y
                    clf.n_estimators -= n_estimators                     
                    hierarchy.set_clf(level,clf,features_used)
                    
                    df = pd.concat(dfs)
                    df.reset_index(inplace=True)
                    df = df.groupby(['index',level+'_gt']).sum()
                    df.reset_index(level=1,inplace=True)
                    train_df = pd.concat(train_dfs)
                    train_df.reset_index(inplace=True)
                    train_df = train_df.groupby(['index']).sum()
                    train_df[level+'_tcounts'] = train_df['training'].astype(int) 
                    train_df[level+'_train'] = train_df['training'].astype(bool)
                    del train_df['training']
                    df[level+'_trains'] = df[level+'_train']
                    del df[level+'_train']
                print(f"Merging data into adata.obsm['{key_added}']")
                total_adata.obsm[key_added] = total_adata.obsm[key_added].merge(df,how='left',left_index=True,
                                          right_index=True,suffixes=['_error',None])   
                total_adata.obsm[key_added] = total_adata.obsm[key_added].merge(train_df,how='left',left_index=True,
                                          right_index=True,suffixes=['_error',None]) 
            if not hierarchy.is_cutoff(level):
                df = pd.DataFrame(index=subset_adata.obs_names)

                print(f'Predicting for {level}...')

                clf,feature_names = hierarchy.get_clf(level)

                if (len(features_used) != len(feature_names)) or (len(features_used) != sum([1 for i,j in zip(feature_names,features_used) if i==j])):
                    logg.debug('Realigning features by feature names...')
                    feature_mask = np.isin(features_used, feature_names)
                    logg.debug(f'sum(feature_mask) = {sum(feature_mask)}, len(feature_names) = {len(feature_names)}')
                    assert sum(feature_mask) == len(feature_names), 'Something is wrong with realignment'
                    subset_X = subset_X[:,feature_mask]
                      
                full_proba, df[level+proba_suffix], df[level+'_class'] = predict_proba_cutoff(subset_X, clf, probability_cutoff)

                if not proba_suffix is None:
                    total_adata.uns[level+proba_suffix] = pd.DataFrame(full_proba,index=np.array(total_adata_obs_names)[df.index.astype(int)],columns=clf.classes_)

            print(f"Merging data into adata.obsm['{key_added}']")
            total_adata.obsm[key_added] = total_adata.obsm[key_added].merge(df,how='left',left_index=True,
                                      right_index=True,suffixes=['_error',None])        
            print('Predicted:')
            print(df[level+'_class'].value_counts())
            if not reference is None:
                logg.debug(f"\n{pd.crosstab(df[level+'_class'],subset_adata.obs[reference]).T.to_string()}")

            if score and (level+'_gt' in total_adata.obsm[key_added].keys()) and (level+'_class' in total_adata.obsm[key_added].keys()):
                plot_confusion(total_adata,level,hierarchy=hierarchy,show=True,key_added=key_added)
                if not batch_key is None:
                    plot_confusion(total_adata,level,hierarchy=hierarchy,show=True,key_added=key_added)
            elif score:
                logg.warning(f'Could not score {level}')
            del subset_X
            del subset_adata
            gc.collect()
        except Exception as e:
            logg.error(f'Ran into an error, skipping {level}',exc_info=1)
        # Attempt to convert to categories and bools to simplify
    cols = total_adata.obsm[key_added].columns
    try:
        cols = total_adata.obsm[key_added].columns
        total_adata.obsm[key_added].loc[:,cols.str.contains('_class')] = \
        total_adata.obsm[key_added].loc[:,cols.str.contains('_class')].astype('category')
        total_adata.obsm[key_added].loc[:,cols.str.contains('_gt')] = \
        total_adata.obsm[key_added].loc[:,cols.str.contains('_gt')].astype('category')
        total_adata.obsm[key_added].loc[:,cols.str.contains('_train')] = \
        total_adata.obsm[key_added].loc[:,cols.str.contains('_train')].astype(bool)
        total_adata.obsm[key_added].loc[:,cols.str.contains('_tcounts')] = \
        total_adata.obsm[key_added].loc[:,cols.str.contains('_tcounts')].astype(float)
    except Exception as e:
        logg.error('Error converting adata.obsm[{key_added}] contents to savable categories',exc_info=1)
    total_adata.obs_names = total_adata_obs_names
    gc.collect()
    return total_adata,hierarchy

def _default_kwargs(kwargs, defaults):
    '''Helper function to combine lists of kwargs, allowing the left (first argument) to override the defaults (second, right argument)'''
    for key in defaults.keys():
        if key not in kwargs.keys():
            kwargs[key] = defaults[key]
    return kwargs

def predict_proba_cutoff(X,clf, probability_cutoff = 0):
    '''Predict using the classifier, using the probability cutoff as described.
    Replace predictions with a ? if they do not match the cutoff —— TODO check that this is actually correctly handled later on.
    X: The matrix of cell features to predict on
    clf: The classifier to predict with
    probability_cutoff: float between 0 and 1, the minimum confidence probability to tolerate before that cells' classification is removed.
    '''
    full_proba = clf.predict_proba(X)
    proba = full_proba
    n_samples = proba.shape[0]
    class_type = clf.classes_.dtype    # all dtypes should be the same, so just take the first
    predictions = np.empty(n_samples, dtype=class_type)

    for k in range(n_samples):
        predictions[k] = clf.classes_[np.argmax(proba[k])]
    
    proba = np.max(proba,axis=1)
    
    for i,prob in enumerate(proba):
        if prob < probability_cutoff:
            predictions[i] = "?"
    return full_proba, proba, predictions
        
def ground_truth(adata,hierarchy,level,data_key,plot=True,reference=None,batch=None,log_reference=True):
    ''' Thresholding for a list of markers, followed by application of these thresholds to ground truths
    Note, defaults to using ADT data, and will only use GEX if ADT data does not exist, or if
    "_gex" is appended to the gene name. If in neither GEX or ADT, checks obs.
    
    TODO fix info
    '''
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
            temp[marker],threshold = thresholding.threshold(marker,adata,data_key,plot=plot,interactive=interactive,
                                                        preset_threshold=preset_threshold,run=True)
            
        except:
            # If thresholding doesn't work, attempt to find it anywhere else
            temp[marker] = utils.get_data(adata,marker,data_key)
            
        if not threshold is None:
            hierarchy.set_threshold(marker,threshold, False, level, batch)
    # Set up and excecute gt clauses
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
                    print(label,clause)
                    print(adata.obs[temp['mask'].values][reference].value_counts())
                del temp['mask']
            clause = "temp.loc["+clause+",'names'] = \'"+label+"\'"
            clauses.append(clause)
    temp['names'] = "?"
    for clause in clauses:
        exec(clause,locals())
    
    templin[level+'_gt'] = temp['names']

    if not reference is None:
        df = pd.crosstab(templin[level+'_gt'],adata.obs[reference]).T
        if log_reference:
            logg.debug('Ground Truth calculated...')
            logg.debug(f"\n{templin[level+'_gt'].value_counts()}")
            logg.debug('\n{}'.format(df.to_string()))
        else:
            print('Ground Truth calculated...')
            print(templin[level+'_gt'].value_counts())
            print(df)
    
    return templin
    
def _get_train_data(X,templin,level,resample_method,max_training,in_danger_noise_checker):
    '''TODO add info'''
    print('Choosing training data...') # Randomly select ground truth dataset for training purposes...80% train, 20% held out
    templin[level+'_train'] = False
    mask = np.random.choice([True,True,True,True,False],size=sum(templin[level+'_gt'] != '?'))
    templin.loc[templin[level+'_gt'] != '?',level+'_train'] = mask
    logg.debug(f"Number of events {len(templin[level+'_train'].to_numpy())}")
    logg.debug(f"X.shape: {X.shape}")
    X = X[templin[level+"_train"].to_numpy()]
    y = templin[templin[level+'_train']][level+'_gt'].values
    for i in set(y):
        # Remove from classes above max_training size
        total = sum(y==i)
        if total > (max_training):
            n_idxs = np.array(range(0,len(y)))
            idxs = n_idxs[y==i]
            idx_remove = np.random.choice(idxs,total-max_training,replace=False)
            idxs_remove = np.array([j in idx_remove for j in n_idxs])
            y = y[~idxs_remove]
            X = X[~idxs_remove]

    logg.info(f"{len(y)} real cells in training set...")
    print('Resampling...')
    X,y = _balance_training_classes(X,y,resample_method,max_training=max_training,in_danger_noise_checker=in_danger_noise_checker)
    logg.debug(f'Selected for training: {len(y)}, with {stats.itemfreq(y)}')
    print(f"Training with {len(y)} events after {resample_method} resampling...")
    assert len(set(y)) > 1, 'Must have more than 1 class to classify...'

    return X,y

def _map_training(X, training_X, obs_names):
    """https://stackoverflow.com/questions/13982804/find-given-row-in-a-scipy-sparse-matrix
        assumes that cell data in X is unique and that this implementation of L2 produces an 
        unique value for a given row of the sparse matrix
    
        Finds number of duplicates of each cell in training data and creates a dataframe of
        the names from X as indices and the number of that cell as the values
        
        X: Anndata.X, sparse matrix of all data with cells as rows and expression level as columns 
        training_X: Anndata.X, subset of X, which is also a sparse matrix, duplicate cells allowed
        obs_names: list like, labels for X object, which are used and indicies of outputted dataframe
        
        Returns:
        Pandas DataFrame, indices are the obs_names and the values are the number of that cell in the 
        training data
        """
    r = range(1,X.shape[1] + 1)
    L2_training_X = training_X.multiply(1/(training_X.multiply(training_X).sum(axis=1))).multiply(r).sum(axis=1)
    L2_X = X.multiply(1/(X.multiply(X).sum(axis=1))).multiply(r).sum(axis=1)
    training_X_to_L2 = Counter(L2_training_X.A.flatten())
    sample_counts = np.vectorize(training_X_to_L2.__getitem__)(L2_X.A.flatten()) 
    return pd.DataFrame.from_dict(dict(zip(obs_names,sample_counts)), orient="index", columns=["training"])


def _in_danger_noise(nn_estimator, samples, y, kind="danger"):
        """Estimate if a set of sample are in danger or noise.
        Adapted from BorderlineSMOTE and SVMSMOTE of the imblearn package
        Parameters
        ----------
        nn_estimator : estimator object
            An estimator that inherits from
            :class:`~sklearn.neighbors.base.KNeighborsMixin` use to determine
            if a sample is in danger/noise.
        samples : {array-like, sparse matrix} of shape (n_samples, n_features)
            The samples to check if either they are in danger or not.
        target_class : int or str
            The target corresponding class being over-sampled.
        y : array-like of shape (n_samples,)
            The true label in order to check the neighbour labels.
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
            
        # Targets can also be labelled "in danger" if they are clustered in a small cluster that's not well-represented
        leiden = nn_estimator.cluster()
        crosstab = pd.crosstab(leiden,y)
        crosstab = ((crosstab.div(crosstab.sum(axis=0),axis=1)) <= .05) & ((crosstab.div(crosstab.sum(axis=1),axis=0)) >= .75) & (crosstab >= 5)
        for i in crosstab:
            for j,k in zip(crosstab[i],crosstab[i].index):
                if j:
                    items[(leiden == k) & (y == i) & (items==0)] = 1
        return items
    
class HVF_PCA_Neighbors():
    def __init__(self,**kwargs):
        # hvf = sklearn.feature_selection.SelectKBest(k=5000)
        # pca = 
        self.neighbors = sklearn.neighbors.NearestNeighbors(**kwargs)
        self.n_neighbors = self.neighbors.n_neighbors
    def fit(self,X):
        '''TODO implement HVG calculation and pca that does not depend on scanpy'''
        adata = anndata.AnnData(X=X,dtype=np.float32)
        sc.pp.highly_variable_genes(adata,n_top_genes = 5000)
        logg.debug('PCA for fit for HVF_PCA_Neighbors')
        sc.tl.pca(adata,n_comps=15)
        self.pca = adata.obsm['X_pca']
        X = self.pca[:,0:15].copy()
        del adata
        self.neighbors = self.neighbors.fit(X)
        return self
    def kneighbors(self,*args,**kwargs):
        return self.neighbors.kneighbors(**kwargs)
    def kneighbors_graph(self,*args,**kwargs):
        return self.neighbors.kneighbors_graph(**kwargs)
    def cluster(self):
        '''TODO implement neighbors and leiden calculations that do not depend on scanpy'''
        adata = anndata.AnnData(X=self.pca,dtype=np.float32)
        adata.obsm['X_pca'] = self.pca
        logg.debug('neighbors for in_danger_noise clustering')
        sc.pp.neighbors(adata,n_pcs=15)
        logg.debug('clustering for in_danger_noise')
        # skip aall this and do adjacency=kneighbors_graph
        sc.tl.leiden(adata,resolution = .5,n_iterations=1)
        if len(adata.obs.leiden.unique()) > min(10,len(adata)/100):
            sc.tl.leiden(adata,resolution = .25,n_iterations=1)
        return adata.obs.leiden.values
    def get_params(self,**kwargs):
        return self.neighbors.get_params(**kwargs)

def _balance_training_classes(X,y, method=None,max_training=None,in_danger_noise_checker=True):
    '''Methods to balance classes, using the imblearn package. Consult their documentation for more details'''
    #methods = ["undersample",'resample','oversample','SMOTE','KMeansSMOTE','BorderlineSMOTE',None]  

    if in_danger_noise_checker or method in ['SMOTE','BorderlineSMOTE','KMeansSMOTE']:
        logg.debug('begining HVF PCA NEIGHBORS')
        k_neighbors = HVF_PCA_Neighbors(n_jobs=-1,n_neighbors=10).fit(X)
        
    if in_danger_noise_checker:
        logg.debug('begining in_danger_noise')
        danger_or_noise = _in_danger_noise(k_neighbors,X,y)
        print(f'Found: {sum(danger_or_noise==2)} noise and {sum(danger_or_noise==1)} in danger of {len(danger_or_noise)} events.')
        print(pd.crosstab(danger_or_noise,y))
        if set(y[danger_or_noise != 2]) == set(y):
            y = y[danger_or_noise != 2]
            X = X[danger_or_noise != 2]
            danger_or_noise = danger_or_noise[danger_or_noise != 2]
        else:
            print('Could not remove noise, as it would remove a category')
        if sum(danger_or_noise) > 0:
            X = sp.vstack([X]+[X[danger_or_noise==1]]*5,format='csr')        
            y = np.hstack([y]+[y[danger_or_noise==1]]*5)
        else:
            print('No "in danger" found')
            
    if method == "oversample":
        sample = imblearn.over_sampling.RandomOverSampler()
    elif method == "resample":
        sample = imblearn.under_sampling.OneSidedSelection()#(sampling_strategy = dic)
    elif method == "undersample":
        sample = imblearn.under_sampling.RandomUnderSampler()
    elif method == "SMOTE":
        sample = imblearn.over_sampling.SMOTE(k_neighbors=k_neighbors,n_jobs=-1)
    elif method == "BorderlineSMOTE":
        sample = imblearn.over_sampling.BorderlineSMOTE(k_neighbors=k_neighbors,n_jobs=-1)
    elif method == "KMeansSMOTE":
        sample = imblearn.over_sampling.KMeansSMOTE(k_neighbors=k_neighbors,n_jobs=-1) 
    if method in ["undersample",'resample','oversample','SMOTE','KMeansSMOTE','BorderlineSMOTE']:
        X, y = sample.fit_resample(X, y)
    if not max_training is None and len(y) > max_training:
        X,y = sklearn.utils.resample(X,y, replace=False,n_samples=max_training)
    return X,y

def terminal_names(adata, obs_column = 'classification',confidence_column = 'conf',conf_drops=None,key_added='lin',voting_reference = None, majority_min=.9):
    '''
    Create a column in the .obs featuring the horizonal concatenation (rightward has priority) of classifications in the key_added
    
    adata: anndata object 
    obs_column: str, column in the .obs to add to the andata containing the classification 
    confidence_column: str, column in the .obs to add to the andata containing the percent of trees in agreement with the classification 
    conf_drops: float, listlike of float, maximum percent of confidence-loss between levels of the hierarchy before an event should be terminated at an earlier subset
    key_added: key in the adata object that contains the classification results
    
    
    voting_reference: str, column in the .obs corresponding to leiden clsuters or overclusters for use in the voting function (if majority voting)
    majority_min: float, between 0 and 1, fraction of events within an individal voting group that must agree for that voted identity to override all memebers of that voting group
    
    returns anndata object with a column for classifications and potentially a column for confidences
    
    '''
    class_names = [x for x in adata.obsm[key_added].columns if "_proba" in x]
    conf_matrix = adata.obsm[key_added].loc[:,class_names]
    class_names = [x for x in adata.obsm[key_added].columns if "_class" in x][1:]
    names_matrix = adata.obsm[key_added].loc[:,class_names].astype(str)
    conf_matrix.columns = names_matrix.columns
    if conf_drops is None:
        conf_drops = -1
    if not isinstance(conf_drops, list):
        conf_drops = [conf_drops]*len(names_matrix.columns)
    confidences = conf_matrix.iloc[:,0]
    confidences.name = None
    labels = names_matrix.iloc[:,0]
    labels.name = None
    blacklisted = conf_matrix.iloc[:,0] < 0
    blacklisted.name = None
    for col, conf_drop in zip(names_matrix.columns,conf_drops):
        mask = conf_matrix.loc[:,col]
        mask.loc[blacklisted] = -1
        mask.name = None
        labels.loc[mask>conf_drop] = names_matrix.loc[mask>conf_drop,col]
        confidences.loc[mask>conf_drop] = conf_matrix.loc[mask>conf_drop,col]
        blacklisted.loc[mask<conf_drop] = True

    adata.obs[obs_column] = labels
    adata.obs[confidence_column] = confidences
    adata.obs[obs_column] = adata.obs[obs_column].astype('category')
    if voting_reference is None:
        return adata
    else:
        def percent_mode(x):
            mode = statistics.mode(x)
            return sum(x == mode)/len(x)
        simple_majority = adata.obs.groupby(voting_reference).agg({obs_column:statistics.mode})[obs_column]
        majority = simple_majority[(adata.obs.groupby(voting_reference).agg({obs_column:percent_mode}) > majority_min)[obs_column]]
        adata.obs['simple_majority_voting'] = adata.obs[voting_reference].map(simple_majority)
        adata.obs['majority_voting'] = adata.obs[voting_reference].map(majority).fillna(adata.obs[obs_column])
        return adata
    
    
    
def idenitfy_group_markers(adata,group1, group2=[], batch_val=None, batch='donor', reference='leiden', key_added='groups', filtered=False, plot=True,
                           min_fold_change=2, min_in_group_fraction=.5, max_out_group_fraction=0.25, n_ups=29, n_downs=30, return_df=True):
    adata.obs[key_added] = 'outgroup'
    group1 = group1 if isinstance(group1,list) else [group1]
    group2 = group2 if isinstance(group2,list) else [group2]
    batch_val = batch_val if isinstance(batch_val,list) else [batch_val]

    if reference == 'leiden':
        group1, group2 = [str(i) for i in group1], [str(i) for i in group2]
        
    group1_mask = adata.obs[reference].isin(group1)
    
    if group2 == []:
        group2_mask = ~adata.obs[reference].isin(group1)
    else:
        group2_mask = adata.obs[reference].isin(group2)
        
    if not batch_val is None and not batch is None:
        batch_mask = adata.obs[batch].isin(batch_val)
        group1_mask, group2_mask = (group1_mask & batch_mask), (group2_mask & batch_mask)
    else:
        batch_mask = [True]*len(adata)
        
    adata.obs.loc[group1_mask,key_added] = 'group1'
    adata.obs.loc[group2_mask,key_added] = 'group2'
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    sc.tl.rank_genes_groups(adata, groupby=key_added, groups=['group1'], reference='group2')
    
    if filtered:
        sc.tl.filter_rank_genes_groups(adata, key='rank_genes_groups',  key_added='rank_genes_groups', groupby='groups',compare_abs=True,
                                       min_fold_change=min_fold_change,min_in_group_fraction=min_in_group_fraction,max_out_group_fraction=max_out_group_fraction)
        
    df=sc.get.rank_genes_groups_df(adata, group='group1')
    df= df[~df.names.isna()]
    if plot:
        if len(df) > n_ups+n_downs:
            ups,downs = df.head(n_ups).names.to_list(), df.tail(n_downs).names.to_list()
        else:
            ups,downs = df.names.to_list(), []
        sc.pl.umap(adata[batch_mask], color=[key_added]+ups+downs, cmap='jet', s=10,ncols=6)
    
    if return_df:
        return df
    return