import gc
import matplotlib.pyplot as plt
import sklearn
import sklearn.calibration
import sklearn.ensemble
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import random
import scanpy as sc
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple, Callable, List
import anndata
import os 

from . import utils
from .logger import logg

def plot_tree(hierarchy, level: str,
              tree_number: int=random.randint(1,10), save: str='tree.png'):
    '''
    Plots a tree from the random forest at a specified level of the classifier. Note these plots are often very unwieldy. See export_graphviz for details.
    
    requires Graphviz which can be obtained with 
    sudo apt-get install graphviz
    
    Parameters
    ----------
    hierarchy
        Hierarchy object which random forest is plotting on
    level
        Name of the level of interest to pull the tree, must have a trained classifier attached to it
    tree_number
        The tree to pick. If undefined, a random tree. 
    save
        The name to save the plot as a .png format. This must be done before displaying the plot as well
    '''
    clf,feature_names = hierarchy.get_clf(level, base=True)
    estimator =  clf.estimators_[tree_number]
    logg.info(clf.estimators_[tree_number].tree_.node_count, 'nodes in tree')
    from sklearn.tree import export_graphviz
    # Export as dot file
    export_graphviz(estimator, out_file='tree.dot', 
                    feature_names = feature_names,
                    class_names = clf.classes_,
                    rounded = True, proportion = False, 
                    precision = 2, filled = True)

    # Convert to png using system command (requires Graphviz)
    from subprocess import call
    call(['dot', '-Tpng', 'tree.dot', '-o', save, '-Gdpi=600'])
    try:
        os.remove('tree.dot')
    except:
        pass
    # Display in jupyter notebook
    import IPython
    plt = IPython.display.Image(filename = save)
    IPython.display.display(plt)
    return

def feature_importances(hierarchy, level: str) -> pd.DataFrame:
    '''Returns a DataFrame of features used in classification and their importances in the random forest of a given level. 
    See sklearn.RandomForestClassifier for more information on interpreting these values.
    
    Parameters
    ----------
    hierarchy
        Hierarchy object containing classification levels
    level
        Name of the level of interest to pull the tree, must have a trained classifier attached to it
        
    Returns
    -------
    df
        Dataframe of features ordered by their importance
    '''
    clf, feature_names = hierarchy.get_clf(level,base=True)
    df = pd.DataFrame(feature_names,columns=['Feature'])
    df['Importance'] = clf.feature_importances_
    df.sort_values('Importance',ascending=False,inplace=True)
    df.reset_index(drop=True,inplace=True)
    return df

def _check_levels(levels: Union[str, List[str]], hierarchy=None, include_cutoff: bool=False) -> List[str]:
    '''
    Checks whether levels provided are valid classification levels of the hierarchy. 
    Otherwise, if hierarchy is defined, the levels provided (a list of str) are filtered by classification levels
    
    Parameters
    ----------
    levels
        Potential levels of the classifier to check the validity of. Levels may be "All" to specify all classification levels if a hierarchy is defined.
    hierarchy
        Object use to check the existence of levels, if not defined all levels are returned
    include_cutoff
        Whether to return all valid levels or just those tagged with is_cutoff
    '''
    if levels == 'All':
        assert not hierarchy is None, 'You must define the hierarchy to use levels of "All"'
        levels = hierarchy.get_classifications()
    elif not isinstance(levels,list):
        levels = [levels]
    
    if hierarchy is None:
        return levels
    else:
        lvls = []
        for level in levels:
                assert level in hierarchy.get_classifications(), f'{level} is not a valid classification layer'
                if include_cutoff or not hierarchy.get_info(level,'is_cutoff'):
                    lvls.append(level)
    return lvls

def _check_pdf_open(save: str) -> Union[str, PdfPages]:
    '''
    For plotting functions that may benefit from saving to a pdf file with pages for each plot
    
    Parameters
    ----------
    save
        Filepath to save pdf to, if pdf filepath, will create and return PdfPages object
    
    Returns
    -------
    pdf_page
        Pdf object that plots to which plots can be saved
    save
        If filepath does not end with .pdf, just returns the save string
    '''
    if not save is None and save.endswith('.pdf'):
        pdf_page = PdfPages(save)
        return pdf_page
    else:
        return save

def _check_pdf_close(save: PdfPages):
    '''
    Check if a pdf object was opened to attempt to close it 
    '''
    if not save is None and not isinstance(save, str):
        save.close()
    return

def _save_or_show(save: Optional[Union[str,PdfPages]], show: bool, dpi: Union[int, str]=300):
    '''
    Save fig and/or show it
    
    Parameters
    ----------
    save
        Object used to save plots
    show
        Whether to display plot after saving it
    dpi
        Resolution to use, see matplotlib.pyplot.savefig for more details
    '''
    if not save is None:
        if isinstance(save,str):
            plt.savefig(save,dpi=dpi)
        else:
            save.savefig(dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close()
    return

def _mask_adata(adata: anndata.AnnData, level: str,
                hc_only: bool=False, holdout_only: bool=False,
                key_added: str='lin') -> anndata.AnnData:
    '''
    Helper function to mask adata functions reliably
    
    Parameters
    ----------
    adata
        Object to apply mask to and from which mask is created
    level
        Level of the classifier to create mask for
    hc_only
        Whether to only include high-confidence data points that were called as pos or neg for the specified level of classification
    holdout_only
        Whether to only include data only non-training data for the specified level of classification. If string supplied,  adata.obsm[key_added][level+holdout_only] is used to mask adata
    key_added
        Key in .obsm where one checks for high-confidence and untrained data 
    
    Returns 
    -------
    adata: AnnData object
        AnnData with any applicable masks potentially applied
    '''
    mask = [True] * adata.n_obs
    if hc_only:
        mask = mask & (adata.obsm[key_added][level+'_hc'] != "?").values & ~(adata.obsm[key_added][level+'_hc']).isna()
    if holdout_only == True:
        mask = mask & (adata.obsm[key_added][level+'_holdout']).values
    elif isinstance(holdout_only, str):
        mask = mask & (adata.obsm[key_added][holdout_only]).values & ~(adata.obsm[key_added][level+'_class']).isna()
    return adata[mask]

def plot_confusion(adata: anndata.AnnData, levels: Union[str, List[str]], 
                   hierarchy=None, key_added: str='lin', 
                   holdout_only: bool=True, batch_key: str=None, 
                   save: str=None, show: bool=True,
                   title_addition: str=None, **kwargs):
    '''
    Determine the performance at a single level by creating a confusion plot using high-confidence thresholds as truth.
    
    Parameters
    ----------
    adata
        Object containing a DataFrame in the .obsm[key_added] specifiying the results of high_confidence thresholding, training, and classification
    levels
        Level or levels of the classifier to create confusion plots for, "all" creates plots for whole hierarchy
    hierarchy
        Hierarchy object with classification level information
    key_added
        Key in .obsm, where information on whether a level had a high-confidence call or was training data
    holdout_only
        Whether to only include data that was not trained on
    batch_key
        Column within the adata.obs that delineates batches
    save
        Filepath to pdf where the plots will be saved
    show
        Whether to show the saved confusion plots
    title_addition
        Phrase to add to the title of the confusion plots
    **kwargs are passed to sklearn.metrics.ConfusionMatrixDisplay.from_predictions()
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, ax = _plot_confusion(adata[batch_mask],level,hierarchy,key_added,holdout_only, **kwargs)
            if not title_addition is None:
                ax.set_title(ax.get_title()+f'\n{title_addition}')
            if not batch_key is None:
                ax.set_title(f'{batch}\n'+ax.get_title())
            _save_or_show(save,show)
    _check_pdf_close(save)
    return

def plot_external_confusion(adata: anndata.AnnData, levels: Union[str, List[str]],
                            hierarchy=None,true_labels: str='_hc', 
                            key_added: str='lin',save: str=None,
                            show: bool=True,title_addition: str=None,
                            **kwargs):
    '''
    Determine the performance by creating a confusion plot using high-confidence thresholds or supplied true labels as truth on data held out from all training processes.
    
    Parameters
    ----------
    adata
        Object containing a DataFrame in the .obsm[key_added] specifiying the results of high_confidence thresholding, training, and classification
    levels
        Level or levels of the classifier to create confusion plots for, "all" creates plots for whole hierarchy
    hierarchy
        Hierarchy object with classification level information
    true_labels
        Column name in adata.obsm[level+true_labels] of "true" event labels, by default high confidence thresholds will be used
    key_added
        Key in .obsm, where information on whether a level had a high-confidence call or was training data
    save
        Filepath to pdf where the plots will be saved
    show
        Whether to show the saved confusion plots
    title_addition
        Phrase to add to the title of the confusion plots
    **kwargs are passed to sklearn.metrics.ConfusionMatrixDisplay.from_predictions()
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    for level in levels:
        fig, ax = _plot_confusion(adata,level,hierarchy,key_added,holdout_only='external_holdout', true_labels=true_labels, **kwargs)
        if not title_addition is None:
            ax.set_title(ax.get_title()+f'\n{title_addition}')
        _save_or_show(save,show)
    _check_pdf_close(save)
    return

def _plot_confusion(adata, level, hierarchy=None, key_added='lin', holdout_only: Union[str,bool]=True, true_labels='_hc', **kwargs):
    '''
    Determine the performance at a single level by creating a confusion plot using high-confidence thresholds or adata.obs[true_labels] as truth.
    
    Parameters
    ----------
    adata
        Object containing high-confidence information and training calls for each event
    level
        Level of the classifier to create a confusion plot for
    hierarchy
        Hierarchy object with classification level information
    key_added
        Key in .obsm, where information on whether a level had a high-confidence call or was training data
    holdout_only
        Whether to only look at the results that were not trained on or if a column is supplied, the column in .obsm[key_added] to check
    true_labels
        name in .obsm[level+true_labels] where "true" cell type labels exist
    **kwargs are passed to sklearn.metrics.ConfusionMatrixDisplay.from_predictions()
    Returns
    -------
    fig, ax: Confusion matrix for the given level of the classifier
    '''
    adata = _mask_adata(adata, level, hc_only=True, holdout_only=holdout_only, key_added=key_added)
    unseen_data = adata.obsm[key_added]
    if true_labels == '_hc':
        unseen_data = unseen_data[(unseen_data[level+true_labels] != '?').values]
    statuses = unseen_data[level+true_labels].unique() 
    (unseen_hc,unseen_cl) = zip(*[(i,j) for i,j in zip(unseen_data[level+true_labels], unseen_data[level+'_class']) if i in statuses]) 
    fig, ax = plt.subplots(figsize=(4,4))
    try:
        labels = hierarchy.subsets_info(level).keys()
        labels = [i for i in labels if i in unseen_hc]
    except:
        labels = sorted(set(unseen_hc))
    try: # Fails on some versions
        sklearn.metrics.ConfusionMatrixDisplay.from_predictions(unseen_hc,unseen_cl,colorbar=False,ax=ax,xticks_rotation='vertical',labels=labels,**kwargs)
    except:
        sklearn.metrics.ConfusionMatrixDisplay(sklearn.metrics.confusion_matrix(unseen_hc,unseen_cl),
                                               colorbar=False,ax=ax,xticks_rotation='vertical',labels=labels,**kwargs)
    p,r,f,s = sklearn.metrics.precision_recall_fscore_support(unseen_hc, unseen_cl,
                                                              average='weighted', zero_division=0)
    ax.set_title(f"{level}, F1 = {f:.2f}\nPrecision = {p:.2f}, Recall = {r:.2f}")
    return fig, ax

def plot_confidence(adata: anndata.AnnData, levels: Union[str, List[str]],
                    hierarchy=None, key_added: str='lin',
                    probability_suffix: str='_probabilities', holdout_only: bool=True,
                    batch_key: str=None, save: str=None,
                    show: bool=True, title_addition: str=None,
                    bins: int=10):
    '''
    Determine how confident classification is for each subset by displaying calibration curves, which compare the events classified at a given class to its confidence.
    
    Parameters
    ----------
    adata
        AnnData object containing .obsm[key_added] and .obs[batch_key]
    levels
        Level or levels of the classifier to create confidence plots for, "all" creates plots for whole hierarchy
    hierarchy
        Hierarchy object with classification level information
    key_added
        Key in .obsm, where results from high-confidence thresholding and predicted probabilities are stored
    probability_suffix
        Suffix in .obsm[key_added] the delineates probability data for the classification
    holdout_only
        Whether to only include data that was not trained on
    batch_key
        Column within the adata.obs that delineates batches
    save
        Filepath to pdf where the plots will be saved
    show
        Whether to show the saved confidence plots
    title_addition
        Phrase to add to the title of the confidence plots
    bins
        Number of bins to split the probability [0,1] interval. See sklearn.calibration.calibration_curve for more details
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, ax = _plot_confidence(adata[batch_mask],level,key_added,probability_suffix,holdout_only,bins=bins)
            if not title_addition is None:
                ax.set_title(ax.get_title()+f'\n{title_addition}')
            if not batch_key is None:
                ax.set_title(f'{batch}\n'+ax.get_title())
            _save_or_show(save,show)
    _check_pdf_close(save)
    return

def _plot_confidence(adata: anndata.AnnData, level: str, 
                     key_added: str='lin', probability_suffix: str='_probability',
                     holdout_only: bool=True, bins: int=10):
    '''
    Plots one vs all confidence thresholds on a specific level of the classifier. 
    Code adapted from: https://stackoverflow.com/questions/58863673/calibration-prediction-for-multi-class-classification
    
    Parameters
    ----------
    adata
        AnnData object to create confidence thresholds on
    level
        Level of the classifier on which to create the thresholds
    key_added
        Key in .obsm, where results from high-confidence thresholding and predicted probabilities are stored
    probability_suffix
        Suffix in .obsm[key_added] the delineates probability data for the classification
    holdout_only
        Whether to only include data that was not trained on
    bins
        Number of bins to split the probability [0,1] interval. See sklearn.calibration.calibration_curve for more details
    Returns
    -------
    fig, ax1: Confidence plot for specified level of the classifier
    '''
    adata = _mask_adata(adata,level,hc_only=True,holdout_only=holdout_only,key_added=key_added)
    df = adata.uns[level+probability_suffix].copy()
    df = df.merge(adata.obsm[key_added][level+'_hc'],how='inner',left_index=True, right_index=True)
    # Plot the Calibration Curve for every class
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")

    for classification in sorted(df[level+'_hc'].unique()):
        if classification == '?' or pd.isna(classification):
            continue
        df[classification+'_OvA'] = (df[level+'_hc'] == classification)
        fraction_of_positives, mean_predicted_value = sklearn.calibration.calibration_curve(df[classification+'_OvA'], df[classification], n_bins=bins,strategy='uniform')
        
        ax1.plot(mean_predicted_value, fraction_of_positives, "s-", label="%s" % (classification, ))
        ax2.hist(df[classification], range=(0, 1), bins=bins, label=classification, histtype="step", lw=2,log=True)
    ax1.set_ylabel("Proportion of samples classified as this class")
    ax1.set_xlabel("Classification confidence")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="lower right")
    ax1.set_title('Calibration plots (reliability curve)')

    ax2.set_xlabel("")
    ax2.set_ylabel("Count (log scale)")
    ax2.legend(loc="upper center", ncol=2)

    plt.tight_layout()
    return fig, ax1

def plot_important_features(adata: anndata.AnnData, levels: Union[str, List[str]],
                            hierarchy, 
                            reference: str='hc -> class', data_key: Optional[Union[str,list]]=utils.DATA_KEY, 
                            holdout_only: bool=False, batch_key: str=None, key_added: str='lin',
                            save: str=None, show: bool=True,
                            title_addition: str=None):
    '''Creates violin plots for the 25 most important genes or proteins for each specified level in levels.
    
    Parameters
    ----------
    adata
        AnnData object to use for plotting
    levels
        Level or levels of the classifier to plot important features of "all" creates plots for whole hierarchy
    hierarchy
        Hierarchy object with classification level information
    reference
        Key in the .obs used to group the violin plots.
        If 'hc -> class', will also add column to obsm containing high-confidence -> level tags to the data
    data_key
        Key in .obsm where feature data exists
    holdout_only
        Whether to only include data that was not trained on
    batch_key
        Column within the adata.obs that delineates batches
    key_added
        Key in .obsm, where results from high-confidence thresholding and predicted probabilities are stored
    save
        Filepath to pdf where the plots will be saved
    show
        Whether to show the saved feature plots
    title_addition
        Phrase to add to the title of the feature plots
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, axs = _plot_important_features(adata[batch_mask],level,hierarchy,reference, holdout_only,data_key)
            axs[0].set_title('GEX')
            if not data_key is None: 
                if isinstance(data_key, list):
                    axs[1].set_title(' + '.join(data_key))
                elif data_key == 'protein':
                    axs[1].set_title('ADT')
                else:
                    axs[1].set_title(data_key)
                if not title_addition is None:
                    axs[1].set_title(axs[1].get_title()+f'{title_addition}\n')
                if not batch_key is None:
                    axs[1].set_title(f'{batch}\n'+axs[1].get_title())
            _save_or_show(save,show)
    _check_pdf_close(save)
    return

def _plot_important_features(adata: anndata.AnnData, level: str,
                             hierarchy, reference: str='hc -> class', holdout_only: bool=False,
                             data_key: Optional[Union[str,list]]=utils.DATA_KEY, key_added: str='lin'):
    '''Creates violin plots for the most important 25 gene, protein, or overall features
    
    Parameters
    ----------
    adata
        AnnData object to use for plotting
    level
        Level of the classifier on which to create the plot of important features 
    hierarchy
        Hierarchy object with classification level information
    reference
        Key in the .obs used to group the violin plots.
        If 'hc -> class', will also add column to obsm containing high-confidence -> level tags to the data
    holdout_only
        Whether to only include data that was not trained on
    data_key
        Key in .obsm where feature data exists
    key_added
        Key in .obsm, where results from high-confidence thresholding and predicted probabilities are stored
    
    Returns
    -------
    fig, ax: Violin plots for the most important 25 gene, protein, or overall features
    '''
    adata = adata.copy()
    if reference == 'hc -> class':
        adata.obs[reference] = adata.obsm[key_added][level+'_hc'].astype(str) +" -> "+ adata.obsm[key_added][level+'_class'].astype(str)
    df = feature_importances(hierarchy, level)
    
    obsm_key = None
    if data_key is not None:
        if not isinstance(data_key, list):
            data_key = [data_key]
        obsm_key = list(set(adata.obsm.keys()).intersection(set(data_key)))
        if len(obsm_key) == 0:
            logg.warn(f'No keys {data_key} in .obsm')
            obsm_key = None
        if len(obsm_key) > 1:
            logg.warn(f'Multiple data_keys in .obsm [{obsm_key}], defaulting to {list(obsm_key)[0]}')
    
    is_protein = df.Feature.str.split('_mod_',expand=True)[1] == obsm_key[0]

    df.Feature = df.Feature.str.split('_mod_',expand=True)[0]
    df_gene = df[~is_protein]
    df_prot = df[is_protein]
    if data_key is None:
        fig, axs = plt.subplots(nrows=1,figsize=(15,60))
        axs = [axs]
    else:
        fig, axs = plt.subplots(nrows=2,figsize=(15,60))

    sc.pl.stacked_violin(adata,df_gene.head(25).Feature.tolist(), reference,swap_axes=True, 
                         dendrogram=True, ax=axs[0], show=False, use_raw=False)
    if data_key is None:
        return fig, axs
    if obsm_key is not None:
        for o_key in obsm_key:
            adata2 = utils.obsm_to_X(adata, data_key=obsm_key)
            sc.pl.stacked_violin(adata2, df_prot.head(25).Feature.tolist(), reference, swap_axes=True,
                         dendrogram=True, ax=axs[1], show=False, use_raw=False)
    else:
        adata2 = utils.obsm_to_X(adata, data_key=obsm_key)
        sc.pl.stacked_violin(adata2, df_prot.head(25).Feature.tolist(), reference, swap_axes=True,
                         dendrogram=True, ax=axs[1], show=False, use_raw=False)
    return fig, axs