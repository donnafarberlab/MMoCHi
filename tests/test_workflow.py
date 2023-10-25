import pytest
import mmochi as mmc
import anndata
import numpy as np 
import pandas as pd 

# To prevent pytests blocking for plot windows
mmc.plotting.plt.show = mmc.plotting.plt.close
mmc.landmarking.plt.show = mmc.landmarking.plt.close
mmc.thresholding.plt.show = mmc.thresholding.plt.close

@pytest.fixture(scope="module")
def files():
    files = ['pbmc_10k_protein_v3.h5','5k_pbmc_protein_v3.h5','pbmc_1k_protein_v3.h5']
    cellranger_versions = ['3.0.0','3.0.2','3.0.0']
    base_url =  f"http://cf.10xgenomics.com/samples/cell-exp/"
    urls = [base_url + f"{v}/{i[:-3]}/{i[:-3]}_filtered_feature_bc_matrix.h5" for i,v in zip(files,cellranger_versions)]
    return (files, urls)

@pytest.fixture(scope="module")
def data_loads(files):
    adatas = mmc.utils.preprocess_adatas(['docs/data/'+file for file in files[0]],backup_urls = files[1],log_CP_ADT=10,log_CP_GEX=100000)
    return adatas
    
def test_exclusive_features(data_loads):
    mmc.utils.generate_exclusive_features(data_loads,'protein')
    mmc.utils.generate_exclusive_features(data_loads,None)
    return 

@pytest.fixture(scope="module")
def data_load(data_loads,files):
    adatas = data_loads
    held_out_batch, held_out_file = adatas.pop(2),files[0].pop(2)
    adata = anndata.concat(adatas,merge='first',keys=files[0],label='batch',index_unique='_')
    return adata

def test_single_data_load_no_url():
    adata = mmc.utils.preprocess_adatas('docs/data/pbmc_10k_protein_v3.h5')
    return adata

def test_single_data_load(files):
    adatas = mmc.utils.preprocess_adatas(['docs/data/'+file for file in files[0]],backup_urls = files[1],log_CP_ADT=10,log_CP_GEX=100000)
    adata = mmc.utils.preprocess_adatas('docs/data/pbmc_10k_protein_v3.h5',backup_urls=files[0],
                                         log_CP_ADT=True,log_CP_GEX=True)
    return adata                                         

def assert_order_preserved(data,new_data,obs,include_zeroes = True):      
    assert len(data) == len(new_data)
    assert all(data.index == new_data.index)
    for i in new_data.columns:
        for j in obs['sample'].unique():
            mask = (data[i] != 0) if not include_zeroes else (data[i] != np.nan)
            data1 = data[(obs['sample']==j) & (mask)]
            df1 = new_data[(obs['sample']==j) & (mask)]
            assert all((sorted(df1.loc[data1.sort_values(i).index][i]) == df1.loc[data1.sort_values(i).index][i])), f'Order of {j} donor not preserved {i}'
    return

def test_save_peak_overrides():
    peak_overrides = {'5k_pbmc_protein_v3.h5':{'CD19':[1.4,2.7]}}
    save_path = 'docs/data/overrides_test.json'
    mmc.save_peak_overrides(save_path,peak_overrides)
    loaded_overrides = mmc.load_peak_overrides(save_path)
    assert peak_overrides == loaded_overrides, 'Loaded overrides are not the same as the saved ones.'
    return save_path

@pytest.mark.parametrize('kwargs', [dict(batch='5k_pbmc_protein_v3.h5', marker='CD19', update_lower=False, update_upper=False),
                                    dict(batch='5k_pbmc_protein_v3.h5', marker='CD19', update_lower=1, update_upper=False),
                                    dict(batch='5k_pbmc_protein_v3.h5', marker='CD19', update_lower=False, update_upper=2),
                                    dict(batch='5k_pbmc_protein_v3.h5', marker='CD19', update_lower=False, update_upper=None),
                                    dict(batch='5k_pbmc_protein_v3.h5', marker='CD19', update_lower=None, update_upper=False),
                                    dict(batch='5k_pbmc_protein_v3.h5', marker='CD19', update_lower=False, update_upper=False,
                                         peak_overrides='docs/data/overrides_test.json'),
                                    dict(batch='cat', marker='CD19', update_lower=2, update_upper=5),
                                    dict(batch='cat', marker='CD19', update_lower=2, update_upper=5, current_peaks=None)
                                    ],
                         ids=['no_change','only_lower','only_upper','single_pos','single_neg','from_save','new_batch','new_batch_no_current'])
def test_update_peak_overrides(kwargs):
    peak_overrides = {'5k_pbmc_protein_v3.h5':{'CD19':[1.4,2.7]}}
    if not 'peak_overrides' in kwargs:
        kwargs['peak_overrides'] = peak_overrides
    if not 'current_peaks' in kwargs:
        kwargs['current_peaks'] = {'5k_pbmc_protein_v3.h5':{'CD19':[1.4,2.7]}}
    test_new_overrides = peak_overrides.copy()
    new_overrides = mmc.update_peak_overrides(**kwargs)
    if kwargs['update_upper'] is False:
        kwargs['update_upper'] = kwargs['current_peaks'][kwargs['batch']][kwargs['marker']][-1]
    if kwargs['update_lower'] is False:
        kwargs['update_lower'] = kwargs['current_peaks'][kwargs['batch']][kwargs['marker']][0]
    
    if not kwargs['batch'] in test_new_overrides:
        test_new_overrides[kwargs['batch']] = dict()
    if kwargs['update_upper'] is None:
        test_new_overrides[kwargs['batch']][kwargs['marker']] = [kwargs['update_lower']]
    else:
        test_new_overrides[kwargs['batch']][kwargs['marker']] = [kwargs['update_lower'],kwargs['update_upper']] 

    assert test_new_overrides == new_overrides, f'Overrides updated wrong {test_new_overrides}, {new_overrides}'
    return



    
@pytest.mark.parametrize("kwargs", [dict(),
                                    dict(show=True),
                                    dict(show=3),
                                    dict(marker_bandwidths={'CD19':0.15,'CD3':0.1}),
                                    dict(single_peaks=['CD19','CD4']),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4,2.7]}}),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[None,2.7]}}),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4]}}),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':['None',1.4]}}),
                                    dict(peak_overrides='docs/data/overrides_test.json')],  
                              ids=['defaults','show','show3','set_marker_bandwidth',
                                   'single_peaks','peak_override','single_positive','single_negative','single_pos_string','peak_overrides_file'])
def test_landmark(data_load,kwargs):
    adata = mmc.landmark_register_adts(data_load,'batch',data_key='protein',**kwargs)
    assert 'landmark_protein' in adata.obsm.keys()
    for i in adata.obsm['protein'].columns:
        assert i in adata.obsm['landmark_protein'].columns
    adata.obs['sample'] = adata.obs.batch
    assert_order_preserved(adata.obsm['protein'], adata.obsm['landmark_protein'],adata.obs, False)

    
@pytest.mark.parametrize("kwargs", [dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4,500]}}),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4,-50]}}) ],ids=["pos", "neg"])
def test_landmark_peak_override_out_range(data_load,kwargs):
    with pytest.raises(IndexError):
        adata = mmc.landmark_register_adts(data_load,'batch',data_key='protein',**kwargs)
    

    
    
@pytest.fixture(scope="module")
def landmarked(data_load):
    adata = mmc.landmark_register_adts(data_load,'batch',data_key='protein')
    return adata

def test_landmark_update(landmarked):
    landmarked.obs['sample'] = landmarked.obs.batch
    peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[None,2.7]}}
    landmarked_new_1 = mmc.update_landmark_register(landmarked.copy(),'5k_pbmc_protein_v3.h5','CD19',peak_overrides,'batch')
    
    landmarked_new_2 = mmc.update_landmark_register(landmarked.copy(),'5k_pbmc_protein_v3.h5','CD19',peak_overrides,'batch',
                                 single_peaks=False,bandwidth={'CD19':0.4})
    assert_order_preserved(landmarked_new_1.obsm['protein'], landmarked_new_1.obsm['landmark_protein'],landmarked_new_1.obs, False)
    assert_order_preserved(landmarked_new_2.obsm['protein'], landmarked_new_2.obsm['landmark_protein'],landmarked_new_2.obs, False)

    assert all(landmarked.obsm['landmark_protein'].loc[landmarked.obs.batch!='5k_pbmc_protein_v3.h5','CD19'] == \
           landmarked_new_1.obsm['landmark_protein'].loc[landmarked_new_1.obs.batch!='5k_pbmc_protein_v3.h5','CD19'])
    assert all(landmarked.obsm['landmark_protein'].loc[landmarked.obs.batch!='5k_pbmc_protein_v3.h5','CD19'] == \
           landmarked_new_2.obsm['landmark_protein'].loc[landmarked_new_2.obs.batch!='5k_pbmc_protein_v3.h5','CD19'])   
    return

def test_landmark_update_zeros(landmarked):
    landmarked = landmarked.copy()
    landmarked.obsm['landmark_protein'].loc[landmarked.obs.batch=='5k_pbmc_protein_v3.h5','CD19'] = 0
    landmark = mmc.update_landmark_register(landmarked.copy(),'5k_pbmc_protein_v3.h5','CD19',{},'batch')
    
    landmarked.obsm['landmark_protein'].loc[landmarked.obs.batch=='5k_pbmc_protein_v3.h5','CD19'] = 0
    landmarked.obsm['landmark_protein'].loc[landmarked.obs.batch=='5k_pbmc_protein_v3.h5','CD19'][0] = 10
    landmark = mmc.update_landmark_register(landmarked.copy(),'5k_pbmc_protein_v3.h5','CD19',{},'batch')
    return

def test_landmark_plots(landmarked):
    mmc.stacked_density_plots(landmarked, landmarked.obsm['protein'].columns, 'batch')
    mmc.stacked_density_plots(landmarked, landmarked.obsm['protein'].columns, None)
    mmc.stacked_density_plots(landmarked, landmarked.obsm['protein'].columns, 'batch',subsample=0.5,save_fig='docs/data/example_plot.pdf')
    mmc.stacked_density_plots(landmarked, 'CD19', 'batch',data_key='protein')
    mmc.stacked_density_plots(landmarked, landmarked.obsm['protein'].columns, 'batch',bw_adjust=3)
    mmc.density_plot(landmarked,'CD19','5k_pbmc_protein_v3.h5','batch','protein')
    mmc.density_plot(landmarked,'CD19','5k_pbmc_protein_v3.h5','batch','landmark_protein')
    mmc.density_plot(landmarked,'CD19','5k_pbmc_protein_v3.h5','batch','protein',bw_adjust = 2)
    mmc.density_plot_total(landmarked,'CD19','5k_pbmc_protein_v3.h5','batch','protein',bw_adjust = 2)

def test_hierarchy_setup():
    h= mmc.Hierarchy(default_max_training = 200)
    h.add_classification('Gross','All','CD3')
    with pytest.raises(AssertionError):
        h.add_subset('Lymphocyte','Gross',dict(neg = ['CD14','CD33','MARCO','MERTK'],any_of=['CD3','CD19','CD127','JCHAIN'],n=1))
    h.add_subset('A','Gross','pos')
    h.add_subset('B','Gross','neg')
    h.add_classification('Mine','A',['CD4','CD8','CD16'])
    h.add_subset('C','Mine',dict(any_of=['CD4','CD8','CD16'],any_ofs_connector='|'))
    h.add_subset('D','Mine',dict(neg=['CD4','CD8','CD16']))
    assert not h.has_clf('Gross')
        
@pytest.fixture(scope="module")
def test_hierarchy():
    h= mmc.Hierarchy(default_min_events=15,default_class_weight = 'balanced', default_max_training = 200)
    h.add_classification('Gross','All', ['CD14','CD33','MARCO','MERTK','CD3','CD19','CD127','JCHAIN']) 
    h.add_subset('Lymphocyte','Gross',dict(neg = ['CD14','CD33','MARCO','MERTK'],any_of=['CD3','CD19','CD127','JCHAIN'],n=1))
    h.add_subset('Myelocyte','Gross',dict(any_of = ['CD14','CD33','MARCO','MERTK'],n=1, neg=['CD3','CD19','CD127','JCHAIN']))
    h.add_classification('Lymphoid','Lymphocyte',['CD3','CD19','CD56','CD127','JCHAIN'])           
    h.add_subset('T cell','Lymphoid',dict(pos=['CD3','CD127'],neg=['CD19','JCHAIN']))
    h.add_subset('B cell','Lymphoid',dict(any_of=['CD19'], neg=['CD3','CD127','JCHAIN']))
    h.add_subset('Plasma cell','Lymphoid',dict(any_of=['JCHAIN'], neg=['CD3','CD127']))
    h.add_subset('NK/ILC','Lymphoid',dict(neg=['CD19','JCHAIN','CD3'], any_of=['CD127','CD56']))
    h.color_dict(True,rot=1,hue=3, mode='DEPTH')
    h.display(True)
    return h

@pytest.fixture(scope="module")
def test_hierarchy2():
    h= mmc.Hierarchy(default_min_events=15,default_class_weight = 'balanced', default_max_training = 200)
    markers = ['CD14','CD33','MARCO','MERTK','CD3','CD19','CD127','JCHAIN'] 
    h.add_classification('Gross','All', markers) 
    h.add_subset('Lymphocyte','Gross',mmc.hc_defs(markers,neg = ['CD14','CD33','MARCO','MERTK'],any_of=['CD3','CD19','CD127','JCHAIN'],n=1))
    h.add_subset('Myelocyte','Gross',mmc.hc_defs(markers,any_of = ['CD14','CD33','MARCO','MERTK'],n=1, neg=['CD3','CD19','CD127','JCHAIN']))
    markers = ['CD3','CD19','CD56','CD127','JCHAIN']
    h.add_classification('Lymphoid','Lymphocyte',markers)           
    h.add_subset('T cell','Lymphoid',mmc.hc_defs(markers,pos=['CD3','CD127'],neg=['CD19','JCHAIN']))
    h.add_subset('B cell','Lymphoid',mmc.hc_defs(markers,any_of=['CD19'], neg=['CD3','CD127','JCHAIN']))
    h.add_subset('Plasma cell','Lymphoid',mmc.hc_defs(markers,any_of=['JCHAIN'], neg=['CD3','CD127']))
    h.add_subset('NK/ILC','Lymphoid',mmc.hc_defs(markers,neg=['CD19','JCHAIN','CD3'], any_of=['CD127','CD56'],any_ofs_connector='|'))
    h.add_classification('batches','T cell',['batch'])
    h.add_subset('batch_1','batches',['one'])
    h.add_subset('batch_2','batches',['two'])
    h.color_dict(True,rot=1,hue=3, mode='DEPTH')
    h.display(True)
    h.display(True, supress_labels = True)
    h.display(False)
    h.get_info()
    return h

def test_threshold(landmarked):
    x = mmc.thresholding.run_threshold('CD19',landmarked,'landmark_protein',(1,2))
    assert type(x) == np.ndarray
    return

# Fails on CLI on newer versions of ipywidgets as there is nowhere to plot
# @pytest.mark.parametrize("mode",['fancy rerun all','fancy fill in'])
# def test_hierarchy_run_thresholds(test_hierarchy,landmarked,mode): 
#     test_hierarchy.run_all_thresholds(landmarked,data_key='landmark_protein',
#                                             batch_key='batch',mode=mode)
#     return test_hierarchy

@pytest.mark.parametrize("mode",['every level','fill in','rerun all'])
def test_hierarchy_run_thresholds2(test_hierarchy,landmarked,mode): 
    test_hierarchy = test_hierarchy.copy()
    test_hierarchy.run_all_thresholds(landmarked,data_key='landmark_protein',
                                            batch_key='batch',mode=mode,limit='CD19',interactive=False,batch_marker_order=True)
    return test_hierarchy

@pytest.fixture(scope="module")
def test_hierarchy_load_thresholds(test_hierarchy):
    hierarchy = test_hierarchy.copy()
    hierarchy.reset_thresholds()
    hierarchy.drop_threshold(None,None,0)
    hierarchy.load_thresholds('docs/data/integrated_thresholds.csv')
    hierarchy.drop_threshold(slice(None),slice(None),0)
    hierarchy.batchless_thresholds()
    hierarchy.load_thresholds('docs/data/integrated_thresholds.csv')
    hierarchy.save_thresholds('docs/data/integrated_thresholds.csv')
    hierarchy.save_thresholds(None)

    assert hierarchy.get_all_markers() == set(['CD14','CD33','MARCO','MERTK','CD3','CD19','CD127','JCHAIN','CD56'])
    return hierarchy

# @pytest.mark.xfail(reason="Issue in sc.pl.embedding cmap handling making it incompatible with matplotlib v3.5.3")
def test_umap_thresh(landmarked, test_hierarchy_load_thresholds):
    print(mmc.__path__)
    import matplotlib
    print(matplotlib.__path__)

    landmarked.obsm['X_umap'] = np.ones((len(landmarked),2))
    mmc.utils.umap_thresh(landmarked,test_hierarchy_load_thresholds,batch_key='batch')   
    
def _verify_lin(adata, classification_levels):
    levels = [i for i in classification_levels if i != 'Removal']
    for level in levels:
        hc_na = adata.obsm['lin'][level+'_hc'].isna()
        holdout_na = adata.obsm['lin'][level+'_holdout'] == False
        tcounts_na = adata.obsm['lin'][level+'_traincounts'].isna()
        train_na = adata.obsm['lin'][level+'_train'] == False
        proba_na = adata.obsm['lin'][level+'_probability'].isna()
        class_na = adata.obsm['lin'][level+'_class'].isna()
        assert all(adata.obsm['lin'][hc_na & holdout_na & tcounts_na & train_na & proba_na & class_na] == adata.obsm['lin'][hc_na])
        assert sum(adata.obsm['lin'][level+'_holdout'][~holdout_na] & adata.obsm['lin'][level+'_train'][~train_na]) == 0
        assert sum((~adata.obsm['lin'][level+'_train'][~holdout_na]) & (adata.obsm['lin'][level+'_traincounts'][~holdout_na] > 0)) == 0
        assert sum((adata.obsm['lin'][level+'_train'][~holdout_na]) & ~(adata.obsm['lin'][level+'_traincounts'][~train_na] > 0)) == 0
        assert sum(["_error" in x for x in adata.obsm['lin'].columns]) == 0, "no failed columns"
        
        
def test_classify_cutoff(landmarked, test_hierarchy_load_thresholds):
    mmc.classifier.DEBUG_ERRORS = True
    test_hierarchy_load_thresholds = test_hierarchy_load_thresholds.copy()
    test_hierarchy_load_thresholds.check_all_markers(landmarked,'landmark_protein')
    test_hierarchy_load_thresholds.default_is_cutoff = True
    mmc.define_external_holdout(landmarked)
    adata,hierarchy = mmc.classify(landmarked, test_hierarchy_load_thresholds.copy(), 'lin', 
                                   'landmark_protein', batch_key='batch',
                                   retrain = True)
    adata = mmc.terminal_names(adata)
    adata = mmc.terminal_names(adata,voting_reference='batch')
    mmc.classifier.DEBUG_ERRORS = False
    return

def test_classify_nparray(landmarked, test_hierarchy_load_thresholds):
    mmc.classifier.DEBUG_ERRORS = True
    landmarked = landmarked.copy()
    landmarked.X = landmarked.X.A
    adata,hierarchy = mmc.classify(landmarked, test_hierarchy_load_thresholds.copy(), 'lin', 
                                   'landmark_protein', batch_key='batch',
                                   retrain = True)
    adata = mmc.terminal_names(adata)
    mmc.classifier.DEBUG_ERRORS = False
    _verify_lin(adata, hierarchy.get_classifications())
    return

def test_classify_defaults(landmarked, test_hierarchy_load_thresholds):
    mmc.classifier.DEBUG_ERRORS = True
    test_hierarchy_load_thresholds.check_all_markers(landmarked,'landmark_protein')
    if 'external_holdout' in landmarked.obsm['lin'].columns:
        landmarked.obsm['lin'].drop('external_holdout',axis=1,inplace=True)
    adata,hierarchy = mmc.classify(landmarked, test_hierarchy_load_thresholds.copy(), 'lin', 
                                   'landmark_protein', batch_key='batch',
                                   retrain = True)
    mmc.define_external_holdout(landmarked)
    adata,hierarchy = mmc.classify(landmarked, test_hierarchy_load_thresholds.copy(), 'lin', 
                                   'landmark_protein', batch_key='batch',
                                   retrain = True)
    adata = mmc.terminal_names(adata)
    adata = mmc.terminal_names(adata,voting_reference='batch')
    hierarchy.save('docs/data/test')
    h = mmc.Hierarchy(load='docs/data/test')
    # assert hierarchy == h, 'Loading does not preserve identity'
    mmc.plot_confusion(adata,'All',hierarchy,show=True,holdout_only=False, save='docs/data/Confusion_plots.pdf')
    mmc.plot_confusion(adata,'All',hierarchy,show=True,holdout_only=True, title_addition='Testing',batch_key='batch')
    mmc.plot_confidence(adata,'All',hierarchy,show=True, save='docs/data/Calibration_plots.pdf')
    mmc.plot_confidence(adata,'All',hierarchy,show=True, title_addition='Testing',batch_key='batch')

    mmc.plot_important_features(adata,'All',hierarchy,batch_key='batch',save='docs/data/ImportantFeatures.pdf')
    mmc.plot_important_features(adata,'All',hierarchy,reference='batch',title_addition='Testing')

    mmc.plot_tree(hierarchy,'Gross',1,save='docs/data/tree.png')
    mmc.classifier.DEBUG_ERRORS = False
    _verify_lin(adata, hierarchy.get_classifications())
    return

def _check_external_holdout(adata, level):
    assert sum(adata.obsm['lin']['external_holdout']) > 0, "nothing marked as external hold out"
    assert sum(adata.obsm['lin']['external_holdout']) != len(adata.obsm['lin']['external_holdout']), "everything marked as external hold out"
    assert sum(adata.obsm['lin']['external_holdout'] & adata.obsm['lin'][level + '_holdout']) == 0, "external and internal hold out overlap"
    assert sum(adata.obsm['lin']['external_holdout'] & adata.obsm['lin'][level + '_opt_holdout']) == 0, "external and opt overlap"
    assert sum(adata.obsm['lin'][level + '_holdout'] & adata.obsm['lin'][level + '_opt_holdout']) == 0, "opt and hold out overlap"
    return
    
def test_various_settings(landmarked,test_hierarchy_load_thresholds):
    mmc.classifier.DEBUG_ERRORS = True
    h = mmc.Hierarchy(load='docs/data/test')
    h.default_min_events=1

    landmarked.var['to_use'] = list([True,True,True,False]*len(landmarked.var_names))[0:len(landmarked.var_names)]
    # print('run 1')
    if 'external_holdout' in landmarked.obsm['lin'].columns:
        landmarked.obsm['lin'].drop('external_holdout',axis=1,inplace=True)
    adata,hierarchy = mmc.classify(landmarked, h.copy(), 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True,features_limit='to_use',
                                   reduce_features_min_cells=0) 
    
    _verify_lin(adata, hierarchy.get_classifications())
    adata.write_h5ad('docs/data/test_adata.h5ad')
    
    h2 = h.copy()
    h2.tree['Gross'].data.features_limit = dict(GEX='All',
                                                landmark_protein=['CD19','CD3','CD20','cat'])
    h2.tree['Lymphoid'].data.features_limit = [i+'_mod_GEX' for i in landmarked.var_names]

    # print('run 2: limit features, external hold out column made but unused')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True)
    _verify_lin(adata, hierarchy.get_classifications())
    h2 = h.copy()
    h2.default_force_spike_ins = ['Lymphocyte']
    # print('run 3: force spike ins')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True,
                                   skip_to='Gross',end_at='Gross')
    _verify_lin(adata, ['Gross',])
    h2 = h.copy()
    h2.default_force_spike_ins = ['Lymphocyte']
    # print('run 3.2: force spike ins external hold out')
    mmc.define_external_holdout(landmarked)
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True,
                                   skip_to='Gross',end_at='Gross', external_holdout=True)
    _verify_lin(adata, ['Gross',])
    _check_external_holdout(adata, 'Gross')
    h2 = h.copy()
    h2.default_calibrate = False
    # print('run 4: no calibration nor external hold out')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True, external_holdout=False)
    _verify_lin(adata, hierarchy.get_classifications())
    # print('run 4.2: no calibration, yes external hold out')
    mmc.define_external_holdout(landmarked)
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True, external_holdout=True)
    _verify_lin(adata, hierarchy.get_classifications())
    assert sum(adata.obsm['lin']['external_holdout']) > 0, "nothing marked as external hold out"
    assert sum(adata.obsm['lin']['external_holdout']) != len(adata.obsm['lin']['external_holdout']), "everything marked as external hold out"
    assert sum(adata.obsm['lin']['external_holdout'] & adata.obsm['lin']['Lymphoid_holdout']) == 0, "external and internal hold out overlap"
    
    h2 = h.copy()
    # print('run 5: no batches')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key=None, retrain = True) 
    _verify_lin(adata, hierarchy.get_classifications())
    mmc.classifier.DEBUG_ERRORS = False
    sh2 = h.copy()
    # print('run 5.2: no batches, external hold out')
    mmc.define_external_holdout(landmarked)
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key=None, retrain = True, external_holdout=True) 
    _verify_lin(adata, hierarchy.get_classifications())
    _check_external_holdout(adata, 'Lymphoid')
    mmc.classifier.DEBUG_ERRORS = False
               
    h2 = h.copy()
    h2.default_min_events = 0.5
    # print('run 6: min events')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', batch_key=None, retrain = True) 
    mmc.classifier.DEBUG_ERRORS = True
    
    h2 = h.copy()
    h2.default_cutoff = True
    # print('run 7: cutoff')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key=None, retrain = True, reference='batch')  
    _verify_lin(adata, hierarchy.get_classifications())
    
    # print('run 8: probability cutoff')
    adata,hierarchy = mmc.classify(landmarked, h.copy(), 'lin', 'landmark_protein', 
                               batch_key='batch', retrain = True, probability_cutoff=0.8) 
    _verify_lin(adata, hierarchy.get_classifications())
               
    h2 = h.copy()
    h2.default_in_danger_noise_checker = 'danger and noise and rebalance'
    # print('run 9: in danger noise and rebalance')
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True)   
    mmc.classifier.DEBUG_ERRORS = False
    h2 = h.copy()
    h2.default_in_danger_noise_checker = 'danger and noise and rebalance'
    # print('run 9.2: in danger noise and rebalance, external hold out')
    mmc.define_external_holdout(landmarked)
    adata,hierarchy = mmc.classify(landmarked, h2, 'lin', 'landmark_protein', 
                                   batch_key='batch', retrain = True, external_holdout=True)   
    _check_external_holdout(adata, 'Lymphoid')
    mmc.classifier.DEBUG_ERRORS = False
    
def test_classify_weights(landmarked, test_hierarchy_load_thresholds):
    test_hierarchy_load_thresholds.check_all_markers(landmarked,'landmark_protein')

    adata,hierarchy = mmc.classify(landmarked, test_hierarchy_load_thresholds.copy(), 'lin', 
                                   'landmark_protein', batch_key='batch',
                                   retrain = True,weight_integration=True)
    adata = mmc.terminal_names(adata)
    adata = mmc.terminal_names(adata,voting_reference='batch')
    return


def test_flatten_children(test_hierarchy):
    test_hierarchy = test_hierarchy.copy()
    test_hierarchy.flatten_children('Lymphocyte')
    test_hierarchy.display()
    return