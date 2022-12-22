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

from . import utils
from .logger import logg # Currently unused

def plot_tree(hierarchy, level: str,
              tree_number: int=random.randint(1,10), save: str='tree.png'):
    '''
    Plots a singular tree from the heirarchy. This can be very useful for demonstrating how the random forest classifier functions
    
    requires Graphviz which can be obtained with 
    sudo apt-get install graphviz
    
    Parameters
    ----------
    hierarchy: Hierarchy class
        Hierarchy object which random forest is plotting on
    level: str
        Name of the level of interest to pull the tree, must have a trained classifier attached to it
    tree_number: int
        The tree to pick. If undefined, a random tree. 
    save: str
        The name to save the plot as a .png format. This must be done before displaying the plot as well
    '''
    clf,feature_names = hierarchy.get_clf(level, base=True)
    estimator =  clf.estimators_[tree_number]
    print(clf.estimators_[tree_number].tree_.node_count, 'nodes in tree')
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
    '''Returns a dataframe of features and their importances in the random forest of a given level. 
    See sklearn.RandomForestClassifier for more information on interpreting these values.
    
    Parameters
    ----------
    hierarchy: Hierarchy class
        Hierarchy object containing classification levels
    level: str
        Name of the level of interest to pull the tree, must have a trained classifier attached to it
        
    Returns
    -------
    df: anndata object
        Dataframe of features ordered by their importance
    '''
    clf, feature_names = hierarchy.get_clf(level,base=True)
    clf.feature_importances_
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
    levels: str or list of str
        Potential levels of the classifier to check the validity of. Levels may be "All" to specify all classification levels. 
    hierarchy: optional, hierarchy object
        Object use to check the existence of levels, if not defined all levels are returned
    include_cutoff: bool
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
    save: str
        Filepath to save pdf to, if pdf filepath, will create and return PdfPages object
    
    Returns
    -------
    pdf_page: PdfPages object
        Pdf object that plots to which plots can be saved
    save: str
        If filepath does not end with .pdf, just returns the save string
    '''
    if not save is None and save.endswith('.pdf'):
        pdf_page = PdfPages(save)
        return pdf_page
    else:
        # logg.warn(f"Filepath must be a string ending in .pdf, returning filepath")
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
    save: PdfPages object
        Object used to save plots
    show: bool
        Whether to display plot after saving it
    dpi: int 
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
                gt_only: bool=False, untrained_only: bool=False,
                key_added: str='lin') -> anndata.AnnData:
    '''
    Helper function to mask adata functions reliably
    
    Parameters
    ----------
    adata: anndata object
        Object to apply mask to and from which mask is created
    level: str
        Level of the classifier to create mask for
    gt_only: bool
        Whether to only include ground truthed data points that were called as pos or neg for the specified level of classification
    untrained_only: bool
        Whether to only include data only non training data for the specified level of classification
    key_added: str
        Key in .obsm where one checks for ground truth and untrained data 
    
    Returns 
    -------
    adata: anndata object
        Object with potentially gt_only and/or untrained_only masks applied
    '''
    mask = [True] * adata.n_obs
    if gt_only:
        mask = mask & (adata.obsm[key_added][level+'_gt'] != "?").values & ~(adata.obsm[key_added][level+'_gt']).isna()
    if untrained_only:
        mask = mask & (adata.obsm[key_added][level+'_train'] != True).values
    return adata[mask]

def plot_confusion(adata: anndata.AnnData, levels: Union[str, List[str]], 
                   hierarchy=None, key_added: str='lin', 
                   hold_out_only: bool=True, batch_key: str=None, 
                   save: str=None, show: bool=True,
                   title_addition: str=None):
    '''
    Determine the performance at a single level by creating confusion plots of performance
    
    Parameters
    ----------
    adata: anndata object
        Object containing ground truthing data and trained calls
    levels: str or list of str
        Level or levels of the classifier to create confusion plots for, "all" creates plots for whole hierarchy
    hierarchy: hierarchy object
        Heirarchy object with classification level information
    key_added: str
        Key in .obsm, where information on whether a level had a ground truth call or was training data
    hold_out_only: bool
        Whether to only include data that was not trained on
    batch_key: str
        Column within the adata.obs that delineates batches
    save: str
        Filepath to pdf where the plots will be saved
    show: bool
        Whether to show the saved confusion plots
    title_addition: str
        Phrase to add to the title of the confusion plots
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, ax = _plot_confusion(adata[batch_mask],level,hierarchy,key_added,hold_out_only)
            if not title_addition is None:
                ax.set_title(ax.get_title()+f'\n{title_addition}')
            if not batch_key is None:
                ax.set_title(f'{batch}\n'+ax.get_title())
            _save_or_show(save,show)
    _check_pdf_close(save)
    return

def _plot_confusion(adata,level,hierarchy=None,key_added='lin',hold_out_only=True):
    '''
    Determine the performance at a single level by creating a confusion plot.
    
    Parameters
    ----------
    adata: anndata object
        Object containing ground truthing information and training calls for each event
    level: str
        Level of the classifier to create a confusion plot for
    hierarchy: hierarchy object
        Heirarchy object with classification level information
    key_added: str
        Key in .obsm, where information on whether a level had a ground truth call or was training data
    hold_out_only: bool
        Whether to only look at the results that were not trained on
    
    Returns
    -------
    fig, ax: Confusion matrix for the given level of the classifier
    '''
    adata = _mask_adata(adata, level, gt_only=True, untrained_only=hold_out_only, key_added=key_added)
    unseen_data = adata.obsm[key_added]
    statuses = unseen_data[level+'_gt'].unique() 
    (unseen_gt,unseen_cl) = zip(*[(i,j) for i,j in zip(unseen_data[level+'_gt'], unseen_data[level+'_class']) if i in statuses]) 
    fig, ax = plt.subplots(figsize=(4,4))
    try:
        labels = hierarchy.subsets_info(level).keys()
        labels = [i for i in labels if i in unseen_gt]
    except:
        labels = sorted(set(unseen_gt))
    try: # Fails on some versions
        sklearn.metrics.ConfusionMatrixDisplay.from_predictions(unseen_gt,unseen_cl,colorbar=False,ax=ax,xticks_rotation='vertical',labels=labels)
    except:
        sklearn.metrics.ConfusionMatrixDisplay(sklearn.metrics.confusion_matrix(unseen_gt,unseen_cl),
                                               colorbar=False,ax=ax,xticks_rotation='vertical',labels=labels)
    p,r,f,s = sklearn.metrics.precision_recall_fscore_support(unseen_gt, unseen_cl,
                                                              average='weighted', zero_division=0)
    ax.set_title(f"{level}, F1 = {f:.2f}\nPrecision = {p:.2f}, Recall = {r:.2f}")
    return fig, ax

def plot_confidence(adata: anndata.AnnData, levels: Union[str, List[str]],
                    hierarchy=None, key_added: str='lin',
                    proba_suffix: str='_proba', hold_out_only: bool=True,
                    batch_key: str=None, save: str=None,
                    show: bool=True, title_addition: str=None,
                    bins: int=10):
    '''
    Performs confidence thresholds on specified levels of the classifier
    
    Parameters
    ----------
    adata: anndata object
    levels: str or list of str
        Level or levels of the classifier to create confidence plots for, "all" creates plots for whole hierarchy
    hierarchy: hierarchy object
        Heirarchy object with classification level information
    key_added: str
        Key in .obsm, where information on ground truth and predicted probabilities exist
    proba_suffix: str
        Suffix in .obsm[key_added] the delineates probability data for the classification
    hold_out_only: bool
        Whether to only include data that was not trained on
    batch_key: str
        Column within the adata.obs that delineates batches
    save: str
        Filepath to pdf where the plots will be saved
    show: bool
        Whether to show the saved confidence plots
    title_addition: str
        Phrase to add to the title of the confidence plots
    bins: int 
        Number of bins to split the probability [0,1] interval. See sklearn.calibration.calibration_curve for more details
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, ax = _plot_confidence(adata[batch_mask],level,key_added,proba_suffix,hold_out_only,bins=bins)
            if not title_addition is None:
                ax.set_title(ax.get_title()+f'\n{title_addition}')
            if not batch_key is None:
                ax.set_title(f'{batch}\n'+ax.get_title())
            _save_or_show(save,show)
    _check_pdf_close(save)
    return

def _plot_confidence(adata: anndata.AnnData, level: str, 
                     key_added: str='lin', proba_suffix: str='_proba',
                     hold_out_only: bool=True, bins: int=10):
    '''
    Performs one vs all confidence thresholds on a specific level of the classifier. 
    Code adapted from: https://stackoverflow.com/questions/58863673/calibration-prediction-for-multi-class-classification
    
    Parameters
    ----------
    adata: anndata object
        anndata object to create confidence thresholds on
    level: str
        Level of the classifier on which to create the thresholds
    key_added: str
        Key in .obsm, where information on ground truth and predicted probabilities exist
    proba_suffix: str
        Suffix in .obsm[key_added] the delineates probability data for the classification
    hold_out_only: bool
        Whether to only include data that was not trained on
    bins: int
        Number of bins to split the probability [0,1] interval. See sklearn.calibration.calibration_curve for more details
    Returns
    -------
    fig, ax1: Confidence plot for specified level of the classifier
    '''
    adata = _mask_adata(adata,level,gt_only=True,untrained_only=hold_out_only,key_added=key_added)
    df = adata.uns[level+proba_suffix].copy()
    df = df.merge(adata.obsm[key_added][level+'_gt'],how='inner',left_index=True, right_index=True)
    # Plot the Calibration Curve for every class
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
    ax2 = plt.subplot2grid((3, 1), (2, 0))
    ax1.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated")

    for classification in sorted(df[level+'_gt'].unique()):
        if classification == '?' or pd.isna(classification):
            continue
        df[classification+'_OvA'] = (df[level+'_gt'] == classification)
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
                            reference: str='gt -> class', data_key: str='landmark_protein', 
                            hold_out_only: bool=False, batch_key: str=None,
                            save: str=None, show: bool=True,
                            title_addition: str=None):
    '''Creates violin plots for the most important 25 gene, protein, or overall features for each specified level in levels
    
    Parameters
    ----------
    adata: anndata object
        anndata object to use for plotting
    levels: str or list of str
        Level or levels of the classifier to plot important features of, "all" creates plots for whole hierarchy
    hierarchy: hierarchy object
        Heirarchy object with classification level information
    reference: str
        Key used to group in violin plots. See scanpy.pl.stacked_violin for more details (groupby param).
        If 'gt -> class', will also add column to obsm containing groundtruth -> level tags to the data
    data_key: str
        Key in .obsm where feature data exists
    hold_out_only: bool
        Whether to only include data that was not trained on
    batch_key: str
        Column within the adata.obs that delineates batches
    save: str
        Filepath to pdf where the plots will be saved
    show: bool
        Whether to show the saved feature plots
    title_addition: str
        Phrase to add to the title of the feature plots
    '''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, axs = _plot_important_features(adata[batch_mask],level,hierarchy,reference, hold_out_only,data_key)
            axs[0].set_title('GEX')
            if not data_key is None:
                if data_key == 'protein':
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
                             hierarchy, reference: str='gt -> class', hold_out_only: bool=False,
                             data_key: str='landmark_protein'):
    '''Creates violin plots for the most important 25 gene, protein, or overall features
    
    Parameters
    ----------
    adata: anndata object
        anndata object to use for plotting
    level: str
        Level of the classifier on which to create the plot of important features 
    hierarchy: hierarchy object
        Heirarchy object with classification level information
    reference: str
        Key used to group in violin plots. See scanpy.pl.stacked_violin for more details (groupby param).
        If 'gt -> class', will also add column to obsm containing groundtruth -> level tags to the data
    hold_out_only: bool
        Whether to only include data that was not trained on
    data_key: str
        Key in .obsm where feature data exists
        
    Returns
    -------
    fig, ax: Violin plots for the most important 25 gene, protein, or overall features
    '''
    adata = adata.copy()
    if reference == 'gt -> class':
        adata.obs[reference] = adata.obsm['lin'][level+'_gt'].astype(str) +" -> "+ adata.obsm['lin'][level+'_class'].astype(str)
    df = feature_importances(hierarchy, level)
    is_protein = df.Feature.str.split('_mod_',expand=True)[1] == data_key
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
    adata2 = utils.obsm_to_X(adata, data_key=data_key)
    sc.pl.stacked_violin(adata2, df_prot.head(25).Feature.tolist(), reference, swap_axes=True,
                         dendrogram=True, ax=axs[1], show=False, use_raw=False)
    return fig, axs
