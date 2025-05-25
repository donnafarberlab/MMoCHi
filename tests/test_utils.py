import pytest
from mmochi import utils
from mmochi import logger
from mmochi.logger import logg
import anndata 
import pandas as pd
import numpy as np

def test_list_tools():
    assert utils._list_tools([1,2,3],'+',[4,5,6]) == [1,2,3,4,5,6]
    assert utils._list_tools(['1','2','3'],'==1') == [True,False,False]
    assert utils._list_tools([1,2,3],'-',[4,5,6,1]) == [2,3]
    assert utils._list_tools([1,2,3],')',[4,5,6,1]) is None

def test_logs():
    logg.info('info')
    logg.debug('debug')
    logg.warn('warn')
    logg.error('error')
    logger._initiate_log()
    logger.log_to_file('docs/data/test_log')
    logg.info('info')
    logg.debug('debug')
    logg.warn('warn')
    logg.error('error')
    return

def test_marker():
    res = utils._marker('cat',['dog','caat','CAT','caTT'],allow_multiple=True)
    assert res == ['CAT','caTT']
    res = utils._marker('cat',['dog','caat','CAT','caTT'],allow_multiple=False)
    assert res == 'CAT'
    res = utils._marker('cat',['dog','caat','CAT_hi','caTT'],allow_multiple=False)
    assert res == 'CAT_hi'
    with pytest.raises(ValueError):
        res = utils._marker('cat',['dog','caat','CATT','caTTT'],allow_multiple=False)

@pytest.fixture(scope="module")
def load_adata():
    adata = utils.preprocess_adatas('docs/data/pbmc_10k_protein_v3.h5',
            backup_urls="http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5")[0]
    return adata

def test_marker_get_data(load_adata):
    # Test all
    adata = load_adata
    utils.marker(adata, 'CD19','protein',False)
    utils.marker(adata, 'CD19',None,False)
    assert utils.marker(adata,['CD19','CD3E','cd4'],None,False) == ['CD19','CD3E','CD4']
    assert set(utils.marker(adata,['CD3D','CD3E'],None,True)) == set(['CD3D','CD3E','CD3EAP'])
    assert all(utils.get_data(adata,'CD19', 'protein') == utils.get_data(adata,'CD19_mod_protein'))
    utils.get_data(adata,'CD19_gex')
    with pytest.raises(ValueError):
        utils.get_data(adata,'CD1999_gex')
    with pytest.raises(AssertionError):
        utils.marker(adata, 'CD1999','protein')
    with pytest.raises(AssertionError):
        utils.marker(adata, 'CD3E','protein')
    with pytest.raises(ValueError):
        utils.get_data(adata, 'CD3E','protein',False,True)
    with pytest.raises(ValueError):
        utils.get_data(adata, 'CD3E','protein',False,True)
    
    utils.get_data(adata,'CD3E','counts',False,True)
    utils.get_data(adata,'CD3e','counts',False,True)
    utils.get_data(adata, 'cd3','protein',False,True)
    assert all(utils.get_data(adata, 'CD3E','protein',False,False) == utils.get_data(adata, 'CD3e_gex','protein',False,False))
    adata.obs = pd.DataFrame('i',index=adata.obs_names,columns=['batch'])
    cols = adata.obs.columns
    assert all(utils.get_data(adata, 'batch_obs') == utils.get_data(adata, 'BATCH_obs'))
    utils.get_data(adata, 'batch_obs',return_source=True)
    with pytest.raises(ValueError):
        utils.get_data(adata, 'CD1999','proteins')  

    