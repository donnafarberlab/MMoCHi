import matplotlib.pyplot as plt
import sklearn
import sklearn.calibration
import sklearn.ensemble
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import random
import scanpy as sc

from . import utils
from .logger import logg # Currently unused

def plot_tree(hierarchy, level,tree_number=random.randint(1,10), plotname='tree.png'):
    '''Plots a singular tree from the heirarchy. This can be very useful for demonstrating how the random forest classifier functions
    hierarchy: Hierarchy class
    level: name of the level of interest to pull the tree, must have a trained classifier attached to it
    tree_number: the tree to pick. If undefined, a random tree. 
    plotname: The name to save the plot as a .png format. This must be done before displaying the plot as well
    requires Graphviz which can be obtained with 
    sudo apt-get install graphviz
    '''
    clf,feature_names = hierarchy.get_clf(level)
    estimator =  clf.estimators_[tree_number]
    print(clf.estimators_[tree_number].tree_.node_count, 'nodes in tree')
    from sklearn.tree import export_graphviz
    # Export as dot file
    export_graphviz(estimator, out_file='tree.dot', 
                    feature_names = features_used,
                    class_names = clf.classes_,
                    rounded = True, proportion = False, 
                    precision = 2, filled = True)

    # Convert to png using system command (requires Graphviz)
    from subprocess import call
    call(['dot', '-Tpng', 'tree.dot', '-o', plotname, '-Gdpi=600'])

    # Display in jupyter notebook
    import IPython
    plt = IPython.display.Image(filename = plotname)
    IPython.display.display(plt)
    return

def feature_importances(hierarchy, level):
    '''Returns a dataframe of features and their importances in the random forest of a given level. 
    See sklearn.RandomForestClassifier for more information on interpreting these values.
    hierarchy: Hierarchy class
    level: name of the level of interest to pull the tree, must have a trained classifier attached to it
    '''
    clf, feature_names = hierarchy.get_clf(level)
    clf.feature_importances_
    df = pd.DataFrame(feature_names,columns=['Feature'])
    df['Importance'] = clf.feature_importances_
    df.sort_values('Importance',ascending=False,inplace=True)
    df.reset_index(drop=True,inplace=True)
    return df

def _check_levels(levels, hierarchy=None,include_cutoff=False):
    '''Checks whether levels provided are valid classification levels of the hierarchy. 
    Levels may be "All" to specify all classification levels. 
    Otherwise, if hierarchy is defined, the levels provided (a list of str) are filtered by classification levels
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
                if include_cutoff or not hierarchy.is_cutoff(level):
                    lvls.append(level)
    return lvls

def _check_pdf_open(save):
    ''' For plotting functions that may benefit from saving to a pdf file with pages for each '''
    if not save is None and save.endswith('.pdf'):
        pdf_page = PdfPages(save)
        return pdf_page
    else:
        return save

def _check_pdf_close(save):
    ''' Check if a pdf was opened to attempt to close it '''
    if not save is None and not isinstance(save, str):
        save.close()
    return

def _save_or_show(save,show,dpi=300):
    ''' Save fig and/or show it'''
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

def _mask_adata(adata, level, gt_only, untrained_only, key_added='lin'):
    '''Helper function to mask adata functions reliably'''
    mask = [True] * adata.n_obs
    if gt_only:
        mask = mask & (adata.obsm[key_added][level+'_gt'] != "?").values & ~(adata.obsm[key_added][level+'_gt']).isna()
    if untrained_only:
        mask = mask & (adata.obsm[key_added][level+'_train'] != True).values
    return adata[mask]

def plot_confusion(adata,levels,hierarchy=None,key_added='lin',hold_out_only=True,batch_key=None,save=None,show=True,title_addition=None):
    '''Determine the performance at a single level'''
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
    '''Determine the performance at a single level'''
    adata = _mask_adata(adata, level, gt_only=True, untrained_only=hold_out_only, key_added=key_added)
    unseen_data = adata.obsm[key_added]
    statuses = unseen_data[level+'_gt'].unique() 
    (unseen_gt,unseen_cl) = zip(*[(i,j) for i,j in zip(unseen_data[level+'_gt'], unseen_data[level+'_class']) if i in statuses]) 
    fig, ax = plt.subplots(figsize=(4,4))
    try: # Fails on some versions
        sklearn.metrics.ConfusionMatrixDisplay.from_predictions(unseen_gt,unseen_cl,colorbar=False,ax=ax,xticks_rotation='vertical')
    except:
        sklearn.metrics.ConfusionMatrixDisplay(sklearn.metrics.confusion_matrix(unseen_gt,unseen_cl),
                                               colorbar=False,ax=ax,xticks_rotation='vertical')
    p,r,f,s = sklearn.metrics.precision_recall_fscore_support(unseen_gt,unseen_cl,average='weighted',zero_division=0)
    ax.set_title(f"{level}, F1 = {f:.2f}\nPrecision = {p:.2f}, Recall = {r:.2f}")
    return fig, ax

def plot_confidence(adata,levels,hierarchy=None,key_added='lin',proba_suffix='_proba',hold_out_only=True, batch_key=None,save=None,show=True,title_addition=None):
    '''TODO'''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, ax = _plot_confidence(adata[batch_mask],level,key_added,proba_suffix,hold_out_only)
            if not title_addition is None:
                ax.set_title(ax.get_title()+f'\n{title_addition}')
            if not batch_key is None:
                ax.set_title(f'{batch}\n'+ax.get_title())
            _save_or_show(save,show)
    _check_pdf_close(save)
    return
def _plot_confidence(adata,level,key_added='lin',proba_suffix='_proba',hold_out_only=True):
    '''Performs one vs all confidence thresholds. 
    Todo add pdf export support & plot vs dont show
    Code adapted from: https://stackoverflow.com/questions/58863673/calibration-prediction-for-multi-class-classification
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
        fraction_of_positives, mean_predicted_value = sklearn.calibration.calibration_curve(df[classification+'_OvA'], df[classification], n_bins=10)

        ax1.plot(mean_predicted_value, fraction_of_positives, "s-", label="%s" % (classification, ))
        ax2.hist(df[classification][df[classification+'_OvA']], range=(0, 1), bins=10, label=classification, histtype="step", lw=2)
    ax1.set_ylabel("Proportion of samples correctly classified as this class")
    ax1.set_xlabel("Average percent of forest in support of the classification")
    ax1.set_ylim([-0.05, 1.05])
    ax1.legend(loc="lower right")
    ax1.set_title('Calibration plots (reliability curve)')

    ax2.set_xlabel("")
    ax2.set_ylabel("Count")
    ax2.legend(loc="upper center", ncol=2)

    plt.tight_layout()
    return fig, ax1

def plot_important_features(adata,levels,hierarchy,key_added='lin',reference='gt -> class', data_key='landmark_protein',hold_out_only=False,batch_key=None,save=None,show=True,title_addition=None):
    '''Right now is looking specific for data_key="protein" and is unreliable in other situations'''
    levels = _check_levels(levels,hierarchy,False)
    save = _check_pdf_open(save)
    batch_masks, batches = utils.batch_iterator(adata,batch_key)
    for batch_mask, batch in zip(batch_masks, batches):
        for level in levels:
            fig, axs = _plot_important_features(adata[batch_mask],level,hierarchy,key_added,reference, hold_out_only,data_key)
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

def _plot_important_features(adata,level,hierarchy,key_added='lin',reference='gt -> class',hold_out_only=False,data_key='landmark_protein'):
    '''Right now is looking specific for data_key="protein" and is unreliable in other situations'''
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

    sc.pl.stacked_violin(adata,df_gene.head(25).Feature.tolist(),reference,swap_axes=True,dendrogram=True,ax=axs[0],show=False,use_raw=False)
    if data_key is None:
        return fig, axs
    adata2 = utils.obsm_to_X(adata, data_key=data_key)
    sc.pl.stacked_violin(adata2,df_prot.head(25).Feature.tolist(),reference,swap_axes=True,dendrogram=True,ax=axs[1],show=False,use_raw=False)
    return fig, axs



### Archived plotting functions ###

# def get_miscalled(adata,level,train='None',key_added='lin'):
#     ''' Returns data that was miscalled between the ground truth and classification steps
#     adata: anndata object
#     level: level of the classifier to probe
#     train: bool, str, None. If True or False, pick data used in training or not, respectively. If str, dont subset by training
#     key_added: location in the .obsm to look for classification information
    
#     Returns anndata objects of miscalled and correctly called data
#     '''
#     if train == True:
#         train = False
#     else:
#         train = True

#     miscalled_data = adata[(adata.obsm[key_added][level+'_gt'] != "?").values &
#                        ~(adata.obsm[key_added][level+'_gt']).isna() &
#                        ~(adata.obsm[key_added][level+'_train']==train) &
#                        (adata.obsm[key_added][level+'_class']!=adata.obsm[key_added][level+'_gt']).values]

#     called_data = adata[(adata.obsm[key_added][level+'_gt'] != "?").values &
#                         ~(adata.obsm[key_added][level+'_gt']).isna() &
#                         ~(adata.obsm[key_added][level+'_train']==train) &
#                         (adata.obsm[key_added][level+'_class']==adata.obsm[key_added][level+'_gt']).values]
#     return miscalled, called

# def miscall_confidence_level(adata,levels,hierarchy=None,key_added='lin',pdf_location = None):
#     ''' Plots histograms of the classifier's confidence of miscalled and correctly called data'''
#     if levels == 'total':
#         levels = hierarchy.get_classifications()
#         print(levels)
#     if not isinstance(levels, list):
#         levels = [levels]    
#     if not pdf_location is None:
#         pp = PdfPages(pdf_location)
#     for level in levels:
#         if hierarchy.is_cutoff(level):
#             print('Skipping cutoff...')
#             continue

#         miscalled_data, called_data = get_miscalled(adata,level,key_added)

#         fig, ax = plt.subplots()
#         called_data.obsm[key_added][level+"_proba"].hist(range=(0,1),bins=40,ax=ax)
#         miscalled_data.obsm[key_added][level+"_proba"].hist(range=(0,1),bins=40,rwidth=.8,ax=ax)
#         ax.set_title(f"Held out {level}: miscall (o) or called (b)")
#         if not pdf_location is None:
#             pp.savefig(fig,dpi=400)
#         plt.show()
#     if not pdf_location is None:
#         pp.close()  
#     return
