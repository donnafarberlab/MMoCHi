import pytest
import mmochi as mmc
import anndata

@pytest.fixture(scope="module")
def data_load():
    files = ['pbmc_10k_protein_v3.h5','5k_pbmc_protein_v3.h5','pbmc_1k_protein_v3.h5']
    cellranger_versions = ['3.0.0','3.0.2','3.0.0']
    base_url =  f"http://cf.10xgenomics.com/samples/cell-exp/"
    urls = [base_url + f"{v}/{i[:-3]}/{i[:-3]}_filtered_feature_bc_matrix.h5" for i,v in zip(files,cellranger_versions)]
    adatas = mmc.utils.preprocess_adatas(['docs/data/'+file for file in files],backup_urls = urls,log_CP_ADT=10,log_CP_GEX=100000)
    held_out_batch, held_out_file = adatas.pop(2),files.pop(2)
    adata = anndata.concat(adatas,merge='first',keys=files,label='batch',index_unique='_')
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
       
@pytest.mark.parametrize("kwargs", [dict(),
                                    dict(show=True),
                                    dict(show=3),
                                    dict(marker_bandwidths={'CD19':0.15,'CD3':0.1}),
                                    dict(single_peaks=['CD19','CD4']),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4,2.7]}})],  
                              ids=['defaults','show','show3','set_marker_bandwidth',
                                   'single_peaks','peak_override'])
def test_landmark(data_load,kwargs):
    adata = mmc.landmark_register_adts(data_load,'batch',data_key='protein',**kwargs)
    assert 'landmark_protein' in adata.obsm.keys()
    for i in adata.obsm['protein'].columns:
        assert i in adata.obsm['landmark_protein'].columns
    adata.obs['sample'] = adata.obs.batch
    assert_order_preserved(adata.obsm['protein'], adata.obsm['landmark_protein'],adata.obs, False)

    
# Test other landmark failures.

@pytest.mark.parametrize("kwargs", [dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4,500]}}),
                                    dict(peak_overrides={'5k_pbmc_protein_v3.h5':{'CD19':[1.4,-50]}}) ],ids=["pos", "neg"])
def test_landmark_peak_override_out_range(data_load,kwargs):
    with pytest.raises(IndexError):
        adata = mmc.landmark_register_adts(data_load,'batch',data_key='protein',**kwargs)
    
    
@pytest.fixture(scope="module")
def landmarked(data_load):
    adata = mmc.landmark_register_adts(data_load,'batch',data_key='protein')
    return adata    

def test_landmark_plots(landmarked):
    mmc.stacked_density_plots(landmarked, landmarked.obsm['protein'].columns, 'batch')
    mmc.stacked_density_plots(landmarked, landmarked.obsm['protein'].columns, 'batch',bw_adjust=3)
    mmc.density_plot(landmarked,'CD19','batch','5k_pbmc_protein_v3.h5','protein')
    mmc.density_plot(landmarked,'CD19','batch','5k_pbmc_protein_v3.h5','protein',bw_adjust = 2)

def test_hierarchy_setup():
    h= mmc.Hierarchy()
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
    h= mmc.Hierarchy(default_min_events=15,default_class_weight = 'balanced')
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

@pytest.mark.parametrize("mode",['fancy rerun all','fancy fill in'])
def test_hierarchy_run_thresholds(test_hierarchy,landmarked,mode):   
    test_hierarchy.run_all_thresholds(landmarked,data_key='landmark_protein',
                                            batch_key='batch',plot_all=False,mode=mode)
    return test_hierarchy

@pytest.fixture(scope="module")
def test_hierarchy_load_thresholds(test_hierarchy):
    hierarchy = test_hierarchy.copy()
    hierarchy.reset_thresholds()
    hierarchy.drop_threshold(None,None,0)
    hierarchy.load_thresholds('docs/data/integrated_thresholds.csv')
    hierarchy.generate_batchless_thresholds
    hierarchy.load_thresholds('docs/data/integrated_thresholds.csv')
    hierarchy.save_thresholds('docs/data/integrated_thresholds.csv')
    hierarchy.save_thresholds(None)

    assert hierarchy.get_all_markers() == set(['CD14','CD33','MARCO','MERTK','CD3','CD19','CD127','JCHAIN','CD56'])
    return hierarchy

def test_classify_defaults(landmarked, test_hierarchy_load_thresholds):
    test_hierarchy_load_thresholds.check_all_markers(landmarked,'landmark_protein')

    adata,hierarchy = mmc.classify(landmarked, test_hierarchy_load_thresholds.copy(), 'lin', 
                                   'landmark_protein', batch_key='batch',
                                   retrain = True)
    adata = mmc.terminal_names(adata)
    adata = mmc.terminal_names(adata,voting_reference='batch')
    hierarchy.save('docs/data/test')
    h = mmc.Hierarchy(load='docs/data/test')
    # assert hierarchy == h, 'Loading does not preserve identity'
    mmc.plot_confusion(adata,'All',hierarchy,show=True,hold_out_only=False, save='docs/data/Confusion_plots.pdf')
    mmc.plot_confidence(adata,'All',hierarchy,show=True, save='docs/data/Calibration_plots.pdf')
    mmc.plot_important_features(adata,'All',hierarchy)
