# Example Hierarchies

A lot of care goes into designing MMoCHi hierarchies, and since their gating strategies are easily distributable. We hope that these hierarchies will provide a good jumping off point, but expect that these will take customization to your datasets (since every dataset will have different cell type compositions and markers may be variably effective across datasets).

We hope to continually update this page with new hierarchies that we or members of the community use! Details on how to submit your hierarchy to this list can be found [below](#submitting-a-mmochi-hierarchy).

## Premade hierarchies

| Name | Short Description | Modalities | Species | Author(s) | Publication |
|:---  |:---               |:---        |:---     |:---       |:---         |
| <a href="#human-t-cell-subsets-v1">Human T cell subsets</a> | Classifier for 8 subsets of αβ T cells and monocytes. | CITE-seq | *Homo sapiens* | Daniel Caron | N/A |
| <a href="#human-immune-subsets-v1">Human immune subsets</a>  | Classifier for over 25 human immune cell subsets across tissue sites. | CITE-seq | *Homo sapiens* | Daniel Caron | N/A |

### Human T cell subsets (v1) 
*Contributed by Daniel Caron*

Classifier for 8 subsets of αβ T cells and monocytes.

CITE-seq: Gene expression (suffixed with `_gex`) and protein (no modifier). Aligned to GRCh38 with Gencode v24 annotation and includes antibodies from a custom universal TotalSeq-A panel of ~270 antibodies (BioLegend: 99786). Classified using all protein-coding genes and all proteins, excluding isotype controls.

*Homo sapiens*: Monocytes and αβ T cells sorted by FACS from human PBMCs.

Additional notes: Uses CD62L in place of CCR7, due to issues with CITE-seq staining of CCR7. CD62L has high concordance with CCR7 in human blood (Sallusto et al, Nature, 1999).
            
```python
h = mmc.Hierarchy()

h.add_classification('Gross','All',['CD14','CD3','CD33']) 
h.add_subset('T cell','Gross',dict(neg = ['CD14','CD33'],pos=['CD3']))
h.add_subset('monocyte','Gross',dict(pos = ['CD14','CD33'],n=1, neg=['CD3']))

h.add_classification('CD4_CD8','T cell',['CD4','CD8'])
h.add_subset('CD4 T cell','CD4_CD8',['pos','neg'])
h.add_subset('CD8 T cell','CD4_CD8',['neg','pos'])

h.add_classification('cd4_mem','CD4 T cell',['CD62L','CD45RA'],clf_kwargs={'max_features':.1})
h.add_subset('cd4_n','cd4_mem',['pos','pos'])
h.add_subset('cd4_cm','cd4_mem',['pos','neg'])
h.add_subset('cd4_em','cd4_mem',['neg','neg',])

h.add_classification('cd8_mem','CD8 T cell',['CD62L','CD45RA'],clf_kwargs={'max_features':.1})
h.add_subset('cd8_n','cd8_mem',['pos','pos'])
h.add_subset('cd8_cm','cd8_mem',['pos','neg'])
h.add_subset('cd8_em','cd8_mem',['neg','neg'])
h.add_subset('cd8_emra','cd8_mem',['neg','pos'])            
```  

### Human immune subsets (v1)
*Contributed by Daniel Caron*

Classifier for over 25 human immune cell subsets across tissue sites.

CITE-seq: Gene expression (suffixed with `_gex`) and protein (no modifier). Aligned to GRCh38 with Gencode v24 annotation and includes antibodies from a custom universal TotalSeq-A panel of ~270 antibodies (BioLegend: 99786). Classified using all protein-coding genes and all proteins, excluding isotype controls.

*Homo sapiens*: CD45+ immune cells magnetically enriched from eight sites across two donors, including lung, airway, lung-associated lymph node, spleen, jejunum epithelial layer, jejunum lamina propria, bone marrow, and blood, using methods optimized for each site (Domínguez Conde, Science, 2022). 

Additional notes: 

- Uses CD62L in place of CCR7, due to issues with CITE-seq staining of CCR7. CD62L has high concordance with CCR7 in human blood (Sallusto et al, Nature, 1999).
- Does not capture non-class-switched memory B cells
- Thresholds for JCHAIN_gex should capture only high expression.

```python
h = mmc.Hierarchy()

h.add_classification('gross','All', ['CD34','KRT7_gex','CD3','CD2','CD19','CD20','JCHAIN_gex','LILRA4_gex','PLD4_gex','CD335','CD33','CD64','OLR1_gex','C1QA_gex',
                                     'HLA-DQA1_gex','HLA-DQB1_gex','CD1c','CD1C_gex','S100A9_gex','S100A8_gex','FCN1_gex','MARCO_gex','MRC1_gex','SPP1_gex','MERTK_gex',
                                     'MPO_gex','ELANE_gex','PRSS57_gex','TPSB2_gex','KIT_gex','CPA3_gex','MZB1_gex','CD352','EGFR_gex','CYTL1_gex','CALD1_gex','COL1A2_gex',
                                     'ESAM_gex','EGFR','Podoplanin','CD326','CD304','MS4A3_gex','KLF1_gex', 'GATA1_gex','TPSAB1_gex', 'MS4A2_gex']) 
h.add_subset('non_immune','gross',dict(any_of=['KRT7_gex','EGFR_gex','CALD1_gex','COL1A2_gex','ESAM_gex','EGFR','Podoplanin','CD326'],n=1,
                                                   any_ofs_connector='|', neg = ['CD3','CD2','CD19','CD20','JCHAIN_gex','CD352','CD335','CPA3_gex','TPSB2_gex','CD304','MS4A3_gex',
                                                                                 'KLF1_gex', 'GATA1_gex','TPSAB1_gex', 'MS4A2_gex']))
h.add_subset('progenitor','gross',dict(any_of=['MPO_gex','ELANE_gex','PRSS57_gex','CYTL1_gex','CD34','KLF1_gex','GATA1_gex','MS4A3_gex'],n=2,
                                                   any_ofs_connector='|', neg = ['CD3','CD2','CD19','CD20','JCHAIN_gex','CD352','CD335','CPA3_gex','TPSB2_gex','KRT7_gex','EGFR_gex',
                                                                                 'CALD1_gex','COL1A2_gex','ESAM_gex','EGFR','Podoplanin','CD326','CD304']))
h.add_subset('lymphocyte','gross',dict(any_of=['CD3','CD2','CD19','CD20','JCHAIN_gex','MZB1_gex','LILRA4_gex','PLD4_gex','CD335','CD352','CD304'], n=1,
                                       neg=['CD34','KRT7_gex','EGFR_gex','CYTL1_gex','TPSB2_gex','CPA3_gex','CD33','CD64','OLR1_gex','C1QA_gex','S100A9_gex','S100A8_gex','MARCO_gex',
                                            'MRC1_gex','SPP1_gex','MPO_gex','ELANE_gex','PRSS57_gex','CALD1_gex','COL1A2_gex','ESAM_gex','EGFR','Podoplanin','CD326','MS4A3_gex',
                                            'KLF1_gex','GATA1_gex','TPSAB1_gex','MS4A2_gex']))

h.add_subset('myelocyte', 'gross', dict(any_of=[['CD33','CD64','OLR1_gex','C1QA_gex'],['HLA-DQA1_gex','HLA-DQB1_gex','CD1c','CD1C_gex'],['S100A9_gex','S100A8_gex','FCN1_gex'],
                                                ['MARCO_gex','MRC1_gex','SPP1_gex','MERTK_gex']], n=[1,3,2,2], 
                                        neg=['CD34','KRT7_gex','CD3','CD2','CD19','CD20','JCHAIN_gex','LILRA4_gex','PLD4_gex','CD335','TPSB2_gex','CPA3_gex','PRSS57_gex','MPO_gex',
                                             'ELANE_gex','EGFR_gex','CYTL1_gex','CALD1_gex','COL1A2_gex','ESAM_gex','EGFR','Podoplanin','CD326','CD304','MS4A3_gex','KLF1_gex', 
                                             'GATA1_gex','TPSAB1_gex','MS4A2_gex'], any_ofs_connector='|'))
h.add_subset('mast_cell', 'gross',dict(neg=['CD34','KRT7_gex','CD3','CD2','CD19','CD20','JCHAIN_gex','LILRA4_gex','PLD4_gex','CD335','MRC1_gex','OLR1_gex','CD64','MPO_gex','ELANE_gex',
                                            'PRSS57_gex','EGFR_gex','CYTL1_gex','CALD1_gex','COL1A2_gex','ESAM_gex','EGFR','Podoplanin','CD326','CD304','MS4A3_gex', 'KLF1_gex'],
                                         pos=['CD33'],any_of=['TPSB2_gex','KIT_gex','CPA3_gex','TPSAB1_gex','MS4A2_gex'],n=1))

h.add_classification('myeloid','myelocyte',['MARCO_gex','MRC1_gex','SEPP1_gex','MERTK_gex','C1QA_gex','S100A9_gex','S100A8_gex','CD64','FCN1_gex','CD1c','CD1C_gex','MPO_gex','ELANE_gex',
                                            'CD14_gex','HLA-DQA1_gex','HLA-DQB1_gex','CD141','CD123'])
h.add_subset('mo_mac','myeloid',dict(any_of=[['MARCO_gex','MRC1_gex','SEPP1_gex','MERTK_gex','C1QA_gex'],['S100A9_gex','S100A8_gex','CD14_gex','CD64','CD123'],['FCN1_gex']],n=[3,3,1],
                                     any_ofs_connector='|', neg=['CD1c','CD1C_gex']))
h.add_subset('dc','myeloid',dict(any_of=['CD1C_gex','CD1c','HLA-DQA1_gex','CD141','HLA-DQB1_gex'],neg=['FCN1_gex','MARCO_gex','MRC1_gex','SEPP1_gex','MERTK_gex','C1QA_gex',
                                                                                               'S100A9_gex','S100A8_gex','CD14_gex','CD64','CD123'],n=2))

h.add_classification('mono_mac','mo_mac',['C1QA_gex','MARCO_gex','MERTK_gex','SEPP1_gex','FCGR3A_gex','S100A9_gex','S100A8_gex','SELL_gex','CD14','CX3CR1','FCN1_gex','MS4A7'])
h.add_subset('macrophage','mono_mac',dict(any_of=['C1QA_gex','MARCO_gex','MERTK_gex','SEPP1_gex'],neg=['FCN1_gex'],pos=['MS4A7'],n=2))
h.add_subset('nc_monocyte','mono_mac',dict(neg=['MARCO_gex','SEPP1_gex','SELL_gex','CD14'],any_of=['FCGR3A_gex','CX3CR1','C1QA_gex'],pos=['FCN1_gex'],n=1))
h.add_subset('c_monocyte','mono_mac',dict(neg=['C1QA_gex','MARCO_gex','MERTK_gex','SEPP1_gex','FCGR3A_gex'], any_of=[['S100A9_gex','S100A8_gex','FCN1_gex'],['SELL_gex','CD14']],n=[2,1]))

h.add_classification('lymphoid','lymphocyte',['CD3','CD19','CD20','MZB1_gex','JCHAIN_gex','CD2_gex','KLRF1_gex','IL7R_gex', 'NCR2_gex','LILRA4_gex','PLD4_gex','TCR_Vd2','TCR_a_b','CD5'],
                     in_danger_noise_checker=False,clf_kwargs={'max_features':0.1})
h.add_subset('t_cell', 'lymphoid', dict(neg = ['CD19','CD20','MZB1_gex','JCHAIN_gex','NCR2_gex','LILRA4_gex','PLD4_gex'], any_of = ['CD3','TCR_a_b','TCR_Vd2','CD5']))
h.add_subset('nk_ilc', 'lymphoid', dict(neg = ['CD3', 'CD19', 'CD20', 'MZB1_gex', 'TCR_Vd2', 'JCHAIN_gex','LILRA4_gex','PLD4_gex','TCR_a_b',"CD5"], 
                                        any_of=['IL7R_gex', 'KLRF1_gex', 'NCR2_gex']))
h.add_subset('b_like', 'lymphoid', dict(any_of = ['CD19', 'CD20', 'MZB1_gex', 'JCHAIN_gex','LILRA4_gex','PLD4_gex'], neg=['CD3', 'IL7R_gex','TCR_Vd2','KLRF1_gex']))

h.add_classification('b_c_like','b_like',['CD19','CD20','JCHAIN_gex','MZB1_gex','LILRA4_gex','PLD4_gex','CD304'])
h.add_subset('b_cell','b_c_like',dict(neg=['MZB1_gex','JCHAIN_gex','LILRA4_gex','CD304'],any_of=['CD19','CD20']))
h.add_subset('plasma','b_c_like',dict(any_of=['MZB1_gex','JCHAIN_gex'],neg=['LILRA4_gex','PLD4_gex','CD20','CD304'],n=2))
h.add_subset('pDC','b_c_like',dict(neg=['CD19','CD20'],any_of=['LILRA4_gex','PLD4_gex','JCHAIN_gex'],pos=['CD304'],n=2))
             
h.add_classification('b_mem','b_cell',['IgD','CD27','IgG'])
h.add_subset('b_naive','b_mem',dict(pos=['IgD'], neg=['CD27','IgG']))
h.add_subset('b_memory','b_mem',dict(any_of=['CD27','IgG']))

h.add_classification('tcr','t_cell',['TRDC_gex','TCR_Vd2','TCR_g_d','TCR_a_b','TRAC_gex'],in_danger_noise_checker=False,clf_kwargs={'max_features':0.1})
h.add_subset('ab_t','tcr',dict(any_of=['TCR_a_b','TRAC_gex'], neg=['TRDC_gex','TCR_Vd2','TCR_g_d']))
h.add_subset('gd_t','tcr',dict(neg=['TCR_a_b'], any_of=['TRDC_gex','TCR_Vd2','TCR_g_d'],n=2))

h.add_classification('cd4cd8','ab_t',['CD4','CD8','CD4_gex','CD8A_gex'],clf_kwargs={'max_features':0.1})
h.add_subset('cd4_t','cd4cd8',dict(pos=['CD4'],neg=['CD8','CD8A_gex']))
h.add_subset('cd8_t','cd4cd8',dict(pos=['CD8','CD8A_gex'],neg=['CD4','CD4_gex']))

h.add_classification('cd4_mem','cd4_t',['FOXP3_gex','CD62L','CD45RA','CCL5_gex','SELL_gex','CD25','CTLA4_gex','CD127','TIGIT_gex','CCR7_gex','CXCR5','PDCD1_gex','BCL6_gex','CXCR3'])
h.add_subset('cd4_naive_cm','cd4_mem',dict(any_of=[['CD62L','CCR7_gex','SELL_gex']],neg=['FOXP3_gex','CCL5_gex','CTLA4_gex'],n=[2], any_ofs_connector='|'))
h.add_subset('cd4_treg','cd4_mem',dict(any_of=[['FOXP3_gex','CTLA4_gex'],['CD25','TIGIT_gex']],neg=['CD127','CD45RA','CD62L','SELL_gex','CCR7_gex'],n=[1,2],any_ofs_connector='|'))
h.add_subset('cd4_effector','cd4_mem',dict(neg=['FOXP3_gex','CD62L','SELL_gex','CCR7_gex','CD45RA','CTLA4_gex'],pos=['CCL5_gex']))

h.add_classification('cd8_mem','cd8_t',['CD62L','CD45RA','CD57','CCL5_gex','SELL_gex','CCR7_gex'])
h.add_subset('cd8_naive_cm','cd8_mem',dict(any_of=['CCR7_gex','CD62L','SELL_gex'], neg=['CD57']))
h.add_subset('cd8_effector','cd8_mem',dict(pos=['CCL5_gex'], neg=['CD62L']))

h.add_classification('cd4_ncm', 'cd4_naive_cm',['CD45RA','CD62L','CD95','CD122','CD45RO','CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'])
h.add_subset('cd4_naive','cd4_ncm',dict(pos=['CD45RA'], any_of=['CD62L'],neg=['CD45RO','CD95','TOX_gex','ICOS_gex','CD122']))
h.add_subset('cd4_cm','cd4_ncm',dict(neg=['CD45RA'],any_of=[['CD62L'],['CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'],['CD95']],n=[1,2,1],any_ofs_connector='|'))

h.add_classification('cd8_eff','cd8_effector',['CD45RA','CCL5_gex'])
h.add_subset('cd8_em','cd8_eff',dict(neg=['CD45RA'],pos=['CCL5_gex']))
h.add_subset('cd8_emra','cd8_eff',dict(pos=['CCL5_gex','CD45RA']))

h.add_classification('cd8_ncm', 'cd8_naive_cm',['CD45RA','CD62L','CD95','CD122','CD45RO','CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'])
h.add_subset('cd8_naive','cd8_ncm',dict(pos=['CD45RA'], any_of=['CD62L'],neg=['CD45RO','CD95','TOX_gex','ICOS_gex','CD122']))
h.add_subset('cd8_cm','cd8_ncm',dict(neg=['CD45RA'],any_of=[['CD62L'],['CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'],['CD95']],n=[1,2,1],any_ofs_connector='|'))

h.add_classification('cd4_res','cd4_effector',['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6'])
h.add_subset('cd4_trm','cd4_res',dict(any_of=['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6'],n=3))
h.add_subset('cd4_tem','cd4_res',dict(neg=['CD69','CD103','CD49a','ITGA1_gex','CXCR6']))
             
h.add_classification('cd8_res','cd8_em',['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6'])
h.add_subset('cd8_trm','cd8_res',dict(any_of=['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6'],n=3))
h.add_subset('cd8_tem','cd8_res',dict(neg=['CD69','CD103','CD49a','ITGA1_gex','CXCR6']))

h.add_classification('nk_ilcs','nk_ilc',['EOMES_gex','GZMH_gex', 'IL7R_gex', 'FCGR3A_gex', 'GZMK_gex', 'GZMB_gex', 'KIT_gex','IL4I1_gex','RORC_gex','NCR2_gex',
                                         'CD335_NKp46','CD103','CD16','CD56'])
h.add_subset('nk_cd56dim','nk_ilcs',dict(any_of=['FCGR3A_gex','GZMB_gex','EOMES_gex','CD16'],n=2, neg=['GZMK_gex','KIT_gex','IL4I1_gex','RORC_gex','NCR2_gex']))
h.add_subset('nk_cd56hi','nk_ilcs',dict(any_of=['GZMK_gex','EOMES_gex','CD56','CD335_NKp46'],n=2, neg=['FCGR3A_gex','KIT_gex','IL4I1_gex','RORC_gex','NCR2_gex','CD103','CD16']))

h.add_subset('ilc_1','nk_ilcs',dict(any_of=['NCR2_gex'], neg=['FCGR3A_gex' ,'EOMES_gex','KIT_gex']))
h.add_subset('ilc_3','nk_ilcs',dict(any_of=['IL7R_gex','KIT_gex','IL4I1_gex','RORC_gex'], neg=['FCGR3A_gex', 'GZMK_gex', 'GZMH_gex', 'GZMB_gex','EOMES_gex','NCR2_gex'],n=2))          
```

## Submitting a MMoCHi hierarchy

If you have designed a hierarchy for your work that you think may be useful for others, we encourage you to submit it for us to share with the community! To do this, open an issue [here](https://github.com/donnafarberlab/MMoCHi/issues/) by clicking the "New issue" button, choose "Submit a hierarchy", and fill out the form!

Note, although pre-trained MMoCHi classifiers can also theoretically be applied across datasets, this is a much more niche application which would require careful handling of marker names across modalities. At this time we are focusing on thresholding strategies, as they are much more easily applied across datasets.