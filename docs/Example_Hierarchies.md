# Example Hierarchies

A lot of care goes into designing MMoCHi hierarchies. Since our hierarchies and gating strategies are easily distributable, we've created this centralized place for sharing them. We hope that these examples will provide a good jumping off point, but expect that these will take customization to apply to your datasets (since every dataset will have different cell type compositions and markers may be variably effective across datasets).

We hope to continually update this page with new hierarchies that we or members of the community use! Details on how to submit your hierarchy to this list can be found [below](#submitting-a-mmochi-hierarchy).

## Tables

<b><a href="#from-mmochi-manuscript">Example Hierarchies from MMoCHi Manuscript</a></b>

These MMoCHi hierarchies were used for analyses in the original MMoCHi manuscript. We have included them here for inspiration and so that users can more easily replicate our analyses.

| Name | Short Description | Modalities | Species | Author(s) | Publication | Training Data | Date Posted
|:---  |:---               |:---        |:---     |:---       |:---         |:---           |:---
| <a href="#human-t-cell-subsets-v1">Human T cell subsets</a> | Classifier for 8 subsets of αβ T cells and monocytes | CITE-seq | *Homo sapiens* | Daniel P. Caron | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) <br/> [[GSE229791]](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229791) | 04JUL23
| <a href="#human-immune-subsets-v1">Human immune subsets</a>  | Classifier for over 25 human immune cell subsets across tissue sites | CITE-seq | *Homo sapiens* | Daniel P. Caron | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) <br/> [[GSE229791]](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229791) | 04JUL23
| <a href="#glioma-v1">Glioma biopsy</a>  | Classifier for neoplastic and non-neoplastic subsets in a high-grade glioma biopsy | scRNA-seq | *Homo sapiens* | Daniel P. Caron | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) | [(Levitin et al., Mol Syst Biol., 2019)](https://doi.org/10.15252/msb.20188557) <br/> [[GSE116621]](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116621) | 23JAN24
| <a href="#sorted-t-and-nk-cells-v1">Ab-seq of sorted T and NK cells</a>  | Classifier for T cell and NK cell subsets adapted for Ab-seq | Ab-seq | *Homo sapiens* | Daniel P. Caron | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) | [(Trzupek et al.,Wellcome Open Res., 2021)](https://doi.org/10.17605/OSF.IO/EDCTN) <br/> [[OSF: EDCTN]](https://doi.org/10.17605/OSF.IO/EDCTN) | 23JAN24
| <a href="#xenium-of-human-lymph-node-v1">10x Xenium of human lymph node</a>  | Classifier for immune and structural cell types in a human lymph node using spatial information | 10x Xenium (Spatial) | *Homo sapiens* | Daniel P. Caron | [(Caron et al., Cell Reports Methods, 2025)](https://doi.org/10.1016/j.crmeth.2024.100938) | [(10x Genomics)](https://www.10xgenomics.com/datasets)<br>[[Human Lymph Node Preview]](https://www.10xgenomics.com/datasets/human-lymph-node-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard) | 23JAN24

<b><a href="#community-submitted">Community Submitted Hierarchies</a></b>

We welcome community members to submit the hierarchies they use to aid others as they design hierarchies for their own datasets. See details [below](#submitting-a-mmochi-hierarchy).

| Name | Short Description | Modalities | Species | Author(s) | Publication | Training Data | Date Posted
|:---  |:---               |:---        |:---     |:---       |:---         |:---           |:---
| <a href="#human-gamma-delta-t-cell-v-delta-subsets-v1">Human γδ T cell Vδ subsets</a> |  Classifier to subset γδ T cells into Vδ1, Vδ2 or 'other' | scRNA-seq | *Homo sapiens* | Joshua I. Gray | [(Gray et al., Science Immunology, 2024)](https://doi.org/10.1126/sciimmunol.adn3954) | [(Gray et al., Science Immunology, 2024)](https://doi.org/10.1126/sciimmunol.adn3954)<br>[[GSE240858]](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240858) | 26JAN24
| <a href="#extended-human-immune-subsets-v1">Extended human immune subsets</a>  | Classifier for 34 human immune subsets across tissue sites | CITE-seq | *Homo sapiens* | Daniel P. Caron & Steven B. Wells | [(Wells et al., bioRxiv, 2024)](https://doi.org/10.1101/2024.01.03.573877) | N/A | 26JAN24

<b><a href="#prior-versions">Prior Versions of Hierarchies</a></b>

Old versions of submitted hierarchies will be saved here.

| Name | Short Description | Modalities | Species | Author(s) | Publication | Training Data | Date Posted
|:---  |:---               |:---        |:---     |:---       |:---         |:---           |:---
| Placeholder  | - | - | - | - | - | - | - |

## Submitting a MMoCHi Hierarchy

If you have designed a hierarchy for your work that you think may be useful for others, we encourage you to submit it for us to share with the community! To do this, open an issue [here](https://github.com/donnafarberlab/MMoCHi/issues/) by clicking the "New issue" button, choose "Submit a hierarchy", and fill out the form! Feel free to share the hierarchies you use regularly, even if they aren't yet published! We can always update links to publications, or post revised versions!

Note, although pre-trained MMoCHi classifiers can also theoretically be applied across datasets, this is a much more niche application which would require careful handling of marker names across modalities. At this time we are focusing on thresholding strategies, as they are much more easily applied across datasets.

## From MMoCHi Manuscript

### Human T cell subsets (v1) 
*Contributed by Daniel P. Caron*

Classifier for 8 subsets of αβ T cells and monocytes.

CITE-seq: Gene expression (suffixed with `_gex`) and protein (no modifier). Aligned to GRCh38 with Gencode v24 annotation and includes antibodies from a custom universal TotalSeq-A panel of ~270 antibodies (BioLegend: 99786). Classified using all protein-coding genes and all proteins, excluding isotype controls.

*Homo sapiens*: Monocytes and αβ T cells sorted by FACS from human PBMCs.

Additional notes: Uses CD62L in place of CCR7, due to issues with CITE-seq staining of CCR7. CD62L has high concordance with CCR7 in human blood (Sallusto et al., Nature, 1999).

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
*Contributed by Daniel P. Caron*

Classifier for over 25 human immune cell subsets across tissue sites.

CITE-seq: Gene expression (suffixed with `_gex`) and protein (no modifier). Aligned to GRCh38 and includes antibodies from a custom universal TotalSeq-A panel of ~270 antibodies (BioLegend: 99786). Classified using all protein-coding genes and all proteins, excluding isotype controls.

*Homo sapiens*: CD45+ immune cells magnetically enriched from eight sites across two donors, including lung, airway, lung-associated lymph node, spleen, jejunum epithelial layer, jejunum lamina propria, bone marrow, and blood, using methods optimized for each site [(Domínguez Conde et al., Science, 2022)]( https://doi.org/10.1126/science.abl5197). 

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

### Glioma (v1)
*Contributed by Daniel P. Caron*

Classifier for neoplastic and non-neoplastic subsets in a high-grade glioma biopsy. 

scRNA-seq paired with averaged chromosome expression: Gene expression (no modifier) and averaged chromosome expression (suffixed with `_aneuploidy`). Aligned to GRCh19. Classified using all protein-coding genes and the ratio of average expression of chromosome 7 to chromosome 10.

*Homo sapiens*: Radiographically guided biopsies were obtained from a patient at the tumor margin and at the tumor core, as described [(Levitin et al., Molecular Systems Biology, 2019)](https://doi.org/10.15252/msb.20188557).

Additional notes: 

- The use of Chr7/Chr10 expression as a tumor-marker is specific to this tumor type and was confirmed in these samples by whole genome sequening. Applying to other cancers will require knowledge of common anueploidy events within those samples.

```python

h = mmc.Hierarchy() 

h.add_classification('Gross','All', ['chr7_chr10_aneuploidy']) 
h.add_subset('Tumor','Gross',dict(pos=['chr7_chr10_aneuploidy']))
h.add_subset('Non-neoplastic','Gross',dict(neg=['chr7_chr10_aneuploidy']))

h.add_classification('Other Populations','Non-neoplastic',['CLDN5','ESAM','TIE1', # Endothelial
                                                           'CD74','C1QB','C1QC','C1QA','HLA-DRA','TYROBP', # Myeloid
                                                           'CNP','PLP1','MBP', # Oligodendrocyte
                                                           'DCN','COL1A2','PDGFRB']) # Pericyte
h.add_subset('Oligodendrocytes','Other Populations',dict(any_of=['CNP','PLP1','MBP'],n=1,
                                                         neg=['CD74','C1QB','C1QC','C1QA','HLA-DRA','TYROBP','CLDN5','ESAM','TIE1','DCN','COL1A2','PDGFRB']))
h.add_subset('Myeloid','Other Populations',dict(any_of=['CD74','C1QB','C1QC','C1QA','HLA-DRA','TYROBP'],n=1,
                                                neg=['CNP','PLP1','MBP','CLDN5','ESAM','TIE1','DCN','COL1A2','PDGFRB']))
h.add_subset('Endothelial','Other Populations',dict(any_of=['CLDN5','ESAM','TIE1'],n=1,
                                                    neg=['CNP','PLP1','MBP','CD74','C1QB','C1QC','C1QA','HLA-DRA','TYROBP','DCN','COL1A2','PDGFRB']))
h.add_subset('Pericytes','Other Populations',dict(any_of=['DCN','COL1A2','PDGFRB'],n=1,
                                                  neg=['CNP','PLP1','MBP','CD74','C1QB','C1QC','C1QA','HLA-DRA','TYROBP','CLDN5','TIE1']))
```

### Sorted T and NK cells (v1)
*Contributed by Daniel P. Caron*

Abridged hierarchy for sorted T and NK cell subsets, derived from the <a href="#human-immune-subsets-v1">Human immune subsets</a> hierarchy, and applied to these Ab-seq data. This was designed to be applied to sorted PBMCs from lupus patients [(Trzupek et al., Wellcome Open Res, 2022)](https://doi.org/10.12688/wellcomeopenres.16883.2).

Ab-seq (BD Biosciences): Gene expression (suffixed with `_gex`) and protein (no modifier), probing for expression of 51 antibodies and a total of 534 unique genes identified to be highly variable across T and NK cell subsets [(See this table)](https://osf.io/fv3x7). Classified using all probed genes and proteins. 

*Homo sapiens*: Six populations of T and NK cells were sorted from PBMCs of patients with SLE. NK: CD56-bright, CD56-dim; T cells: CD8+, CD4+ CD127hi CD25- (CD127hi Tconv),  CD4+ CD127lo CD25- (CD127lo Tconv), CD4+ CD127lo CD25++ (CD4+ Tregs) [(Trzupek et al., Wellcome Open Res, 2022)](https://doi.org/10.12688/wellcomeopenres.16883.2).

Additional notes: 

- In these data, proliferating cells are classified into their respective populations, but could be further delineated.

```python

h = mmc.Hierarchy() 

h.add_classification('lymphoid','All',['CD3','CD5','IL7R_gex','KLRF1','CD56','CD127'],
                     in_danger_noise_checker=False,clf_kwargs={'max_features':0.1})
h.add_subset('t_cell', 'lymphoid', dict(any_of = ['CD3','CD5']))
h.add_subset('nk', 'lymphoid', dict(neg = ['CD3',"CD5"], any_of=['IL7R_gex', 'KLRF1','CD56','CD127']))

h.add_classification('cd4cd8','t_cell', ['CD4','CD8','CD4_gex','CD8A_gex'],clf_kwargs={'max_features':0.1})
h.add_subset('cd4_t','cd4cd8', dict(pos=['CD4'],neg=['CD8','CD8A_gex']))
h.add_subset('cd8_t','cd4cd8', dict(pos=['CD8','CD8A_gex'],neg=['CD4','CD4_gex']))

h.add_classification('cd4_mem','cd4_t', ['FOXP3_gex','CD62L','CD45RA','CCL5_gex','SELL_gex','CD25','CTLA4_gex','CD127','TIGIT_gex','CCR7_gex'])
h.add_subset('cd4_naive_cm','cd4_mem', dict(any_of=[['CD62L','CCR7_gex','SELL_gex']],neg=['FOXP3_gex','CCL5_gex','CTLA4_gex'],n=[2], any_ofs_connector='|'))
h.add_subset('cd4_treg','cd4_mem', dict(any_of=['FOXP3_gex','CTLA4_gex','CD25','TIGIT_gex'],neg=['CD127','CD45RA'],n=1))
h.add_subset('cd4_effector','cd4_mem', dict(neg=['FOXP3_gex','CD62L','SELL_gex','CCR7_gex','CD45RA','CTLA4_gex'],pos=['CCL5_gex']))

h.add_classification('cd8_mem','cd8_t', ['CD62L','CD45RA','CCL5_gex','SELL_gex','CCR7_gex'])
h.add_subset('cd8_naive_cm','cd8_mem', dict(any_of=['CCR7_gex','CD62L','SELL_gex']))
h.add_subset('cd8_effector','cd8_mem', dict(pos=['CCL5_gex'], neg=['CD62L']))

h.add_classification('cd4_ncm','cd4_naive_cm', ['CD45RA','CD62L','CD95','CD45RO','CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'])
h.add_subset('cd4_naive','cd4_ncm', dict(pos=['CD45RA'], any_of=['CD62L'],neg=['CD45RO','CD95','TOX_gex','ICOS_gex']))
h.add_subset('cd4_cm','cd4_ncm', dict(neg=['CD45RA'],any_of=[['CD62L'],['CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'],['CD95']],
                                     n=[1,2,1],any_ofs_connector='|'))

h.add_classification('cd8_eff','cd8_effector', ['CD45RA','CCL5_gex'])
h.add_subset('cd8_em','cd8_eff', dict(neg=['CD45RA'],pos=['CCL5_gex']))
h.add_subset('cd8_emra','cd8_eff', dict(pos=['CCL5_gex','CD45RA']))

h.add_classification('cd8_ncm', 'cd8_naive_cm', ['CD45RA','CD62L','CD95','CD45RO','CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'])
h.add_subset('cd8_naive','cd8_ncm', dict(pos=['CD45RA'], any_of=['CD62L'],neg=['CD45RO','CD95','TOX_gex','ICOS_gex']))
h.add_subset('cd8_cm','cd8_ncm', dict(neg=['CD45RA'],any_of=[['CD62L'],['CXCR5','PDCD1_gex','BCL6_gex','CXCR3','TOX_gex','ICOS_gex'],['CD95']],
                                     n=[1,2,1],any_ofs_connector='|'))

h.add_classification('nks','nk', ['EOMES_gex','GZMH_gex', 'IL7R_gex', 'FCGR3A_gex', 'GZMK_gex', 'GZMB_gex', 'KIT_gex','CD335','CD16','CD56'])
h.add_subset('nk_cd56dim','nks', dict(any_of=['FCGR3A_gex','GZMB_gex','EOMES_gex','CD16'],n=2, neg=['GZMK_gex','KIT_gex']))
h.add_subset('nk_cd56hi','nks', dict(any_of=['GZMK_gex','EOMES_gex','CD56','CD335'],n=2, neg=['FCGR3A_gex','KIT_gex','CD16']))

```
### Xenium of human lymph node (v1)
*Contributed by Daniel P. Caron*

Classifier for immune and structural cell types in a human lymph node using spatial information.

10x Xenium data: paired gene expression by probing (no modifier) and morphological features extracted from the imaging data (suffixed with `_physical_attributes`). Classified using all probed genes (with controls and unused probes removed) and calculated cellular and nuclear area and circularity, as well as nuclear DAPI intensity (see our manuscript for more details). Gene expression was probed using the 377 genes in the [Xenium Human Multi-Tissue and Cancer Panel](https://www.10xgenomics.com/support/in-situ-gene-expression/documentation/steps/panel-design/pre-designed-xenium-gene-expression-panels)

*Homo sapiens*: Adult human non-diseased lymph node, processed from FFPE tissue using the 10x protocols described [here](https://www.10xgenomics.com/datasets/human-lymph-node-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard). 

Additional notes: 

- The proliferating B cells profiled here are likely mostly germinal center B cells, but not all the markers necessary (e.g. *MME*) are present for their robust identification
- High amounts of background expression on markers such as *CD4* and *CD3E* require slightly higher thresholds

```python
h = mmc.Hierarchy()

h.add_classification('Gross','All', ['PTPRC','CD163','CD14','CD3E','CD19','MS4A1','EGFL7','CD34','KRT7','MARCO','LILRA4','GZMB','PDGFRA', 'PDGFRB','PLD4','CLEC10A','TRAC','MZB1']) 
h.add_subset('Lymphocyte','Gross',dict(neg = ['CD163','EGFL7','KRT7','MARCO','CD14','CLEC10A'],any_of=['CD3E','CD19','MS4A1','LILRA4','GZMB','PLD4','TRAC','PTPRC','MZB1'],n=2))
h.add_subset('Myeloid','Gross',dict(any_of = ['CD163','CD14','MARCO','CLEC10A'],n=1, neg=['CD3E','CD19','MS4A1','EGFL7','KRT7','LILRA4','GZMB','MZB1']))
h.add_subset('Non-Immune','Gross',dict(neg = ['PTPRC','CD14','CD3E','CD163','CD19','MS4A1','MARCO','LILRA4','GZMB','CLEC10A','MZB1'],n=[1,2], any_of = [['EGFL7','CD34','KRT7'],['PDGFRA','PDGFRB']],
                                       any_ofs_connector='|'))

h.add_classification('Structural','Non-Immune', ['EGFL7','PDGFRA','PDGFRB','PCOLCE','VWF']) 
h.add_subset('Stromal cell','Structural', dict(neg = ['VWF','EGFL7'],any_of=[['PDGFRA','PDGFRB'],['PCOLCE']],any_ofs_connector='|',n=[2,1]))
h.add_subset('Endothelial cell','Structural', dict(any_of = ['VWF','EGFL7'],neg=['PCOLCE','PDGFRB']))

h.add_classification('Mac/DC','Myeloid', ['CD1C','CLEC10A','MS4A4A','MARCO','CD14','CD163']) 
h.add_subset('Macrophage','Mac/DC', dict(neg = ['CD1C','CLEC10A'],any_of=['MS4A4A','MARCO','CD14','CD163']))
h.add_subset('DC','Mac/DC', dict(any_of = ['CD1C','CLEC10A'],neg=['MS4A4A','MARCO','CD14','CD163']))

h.add_classification('Lymphoid','Lymphocyte', ['CD3E','CD3D','CD247','TRAC','CD19','MS4A1','MZB1','LILRA4','GZMB','IRF8', 'PLD4'],clf_kwargs=dict(max_features=0.1))           
h.add_subset('T cell','Lymphoid', dict(any_of=['CD3E','CD3D','CD247','TRAC'],n=2,neg=['CD19','MS4A1','MZB1','LILRA4','GZMB','PLD4']))
h.add_subset('B cell','Lymphoid', dict(any_of=['CD19','MS4A1','MZB1'], neg=['CD3E','CD3D','CD247','TRAC','LILRA4','GZMB']))
h.add_subset('pDC','Lymphoid', dict(any_of=['LILRA4','GZMB','IRF8'], neg=['CD3E','CD3D','CD247','TRAC','CD19','MS4A1']))

h.add_classification('B_mem','B cell', ['MKI67','PCNA','CD27','TNFRSF13B','TCL1A','MZB1','PRDM1','FAS','CD19','MS4A1'])
h.add_subset('Naive B cell','B_mem', dict(neg=['CD27','MKI67','PCNA'],any_of=['TCL1A']))
h.add_subset('Plasma cell','B_mem', dict(neg=['TCL1A','MKI67','PCNA'],any_of=['MZB1','PRDM1']))

h.add_subset('Memory B cell','B_mem', dict(any_of=['TNFRSF13B','CD27','FAS'],neg=['TCL1A','MKI67','PCNA']))
h.add_subset('Plasmablast','B_mem', dict(any_of=[['MKI67','PCNA'],['MZB1','PRDM1']],neg=['CD19'],n=[1,1]))
h.add_subset('Prolif. B cell','B_mem', dict(any_of=[['MKI67','PCNA'],['CD19','MS4A1']],neg=['MZB1','PRDM1'],n=[1,1]))

h.add_classification('CD4_CD8','T cell', ['CD4','FOXP3','CD8A','GZMB'])
h.add_subset('CD4 T cell','CD4_CD8', dict(any_of=['CD4','FOXP3'],neg=['CD8A','GZMB']))
h.add_subset('CD8 T cell','CD4_CD8', dict(neg=['CD4','FOXP3'],any_of=['CD8A','GZMB']))

h.add_classification('CD4_subsets','CD4 T cell', ['IL2RA','FOXP3','CCL5','SELL','MKI67','PCNA','CTLA4'])
h.add_subset('CD4 Treg','CD4_subsets', dict(any_of=['IL2RA','FOXP3','CTLA4'],n=2,neg=['MKI67','PCNA']))
h.add_subset('CD4 N/CM','CD4_subsets', dict(neg=['IL2RA','FOXP3','MKI67','PCNA','CCL5'],pos=['SELL']))
h.add_subset('CD4 EM','CD4_subsets', dict(neg=['IL2RA','FOXP3','MKI67','PCNA','SELL'],pos=['CCL5']))
h.add_subset('Prolif. CD4 T cell','CD4_subsets', dict(any_of=['MKI67','PCNA']))

h.add_classification('CD8_subsets','CD8 T cell', ['CCL5','SELL','MKI67','PCNA'])
h.add_subset('CD8 N/CM','CD8_subsets', dict(neg=['MKI67','PCNA','CCL5'],pos=['SELL']))
h.add_subset('CD8 EM','CD8_subsets', dict(neg=['MKI67','PCNA','SELL'],pos=['CCL5']))
h.add_subset('Prolif. CD8 T cell','CD8_subsets', dict(any_of=['MKI67','PCNA']))
```

## Community Submitted

### Human gamma delta T cell V delta Subsets (v1) 
*Contributed by Joshua I. Gray*

Classifer for the subsetting γδ T cells into Vδ1, Vδ2 or 'other'
 
scRNA-seq: Gene expression (no modifier). Aligned to GRCh38 with Gencode v24 annotation

*Homo sapiens*:  Pediatric γδ T cells sorted from spleen, lung, lung-associated lymph node, jejunum lamina propria and mesenteric lymph nodes. Additional adult samples from equivalent tissues extracted from [(Caron et al., bioRxiv, 2023)](https://doi.org/10.1101/2023.07.06.547944) <br/> [[GSE229791]](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229791)

Additional notes:  With additional cells this could be fine tuned to include less common Vδ subsets (e.g.) Vδ3, as well as adding Vγ chain classification to fully elucidate the TCR pairings of γδ T cell subsets. Addition of scTCR-seq or CITE-seq data would further improve this resolution.

```python
h = mmc.Hierarchy()

h.add_classification('tcr','All',['TRDV1','TRDV2','TRDV3','TRAV38-2DV8','TRAV23DV6','TRAV29DV5','TRAV36DV7','TRAV14DV4'],)
h.add_subset('vd1','tcr',dict(pos=['TRDV1'],neg=['TRDV2','TRDV3', 'TRAV38-2DV8','TRAV23DV6','TRAV29DV5','TRAV36DV7','TRAV14DV4']))
h.add_subset('vd2','tcr',dict(pos=['TRDV2'],neg=['TRDV1','TRDV3','TRAV38-2DV8','TRAV23DV6','TRAV29DV5','TRAV36DV7','TRAV14DV4']))
h.add_subset('other','tcr',dict(any_of=['TRDV3','TRAV38-2DV8','TRAV23DV6','TRAV29DV5','TRAV36DV7','TRAV14DV4'],neg=['TRDV1','TRDV2']))
```  

### Extended human immune subsets (v1)
*Contributed by Daniel P. Caron and Steven B. Wells*

Classifier for 34 human immune cell subsets across tissue sites. 

CITE-seq: Gene expression (suffixed with `_gex`) and protein (no modifier). Aligned to GRCh38 with Gencode v24 annotation and includes antibodies from the universal TotalSeq-C panel of ~130 antibodies (BioLegend: 399905). Classified using all protein-coding genes and all proteins, excluding isotype controls.

*Homo sapiens*: CD45+ immune cells magnetically enriched from >10 sites across 24 donors, including lung, airway, lung-associated lymph node, spleen, jejunum epithelial layer, jejunum lamina propria, bone marrow, and blood, using methods optimized for each site (Wells et al., bioRxiv, 2024). 

Additional notes: 

- All non-immune cells MUST be removed from the dataset prior to classification. These subsets can often be easily identified using manual cluster-based annotation. 
- Thresholds for JCHAIN_hi, CD33_hi, HLA-DQA1_gex, and HLA-DQB1_gex should capture only high expression. 
- In danger noise checking was modified/disabled based on whether the cell types clustered separately.
- Assumes TRM share a similar expression profile across tissues, but may be more precise for TRM identification if trained on individual tissues or if TRM from each tissue are separated.

```python

h = mmc.Hierarchy()

h.add_classification('gross','All', ['CD2_gex','CD3','CD127','CD19','CD20','CD163','CD33','CD64','CLEC9A_gex','ELANE_gex','JCHAIN_gex',
                                     'LAMP3_gex','MPO_gex','MRC1_gex','OLR1_gex','CLEC12A_gex','PRSS57_gex','CYTL1_gex','SELENOP_gex','TPSB2_gex','CCL19_gex',
                                     'LILRA4_gex','PLD4_gex','MMP9_gex','G0S2_gex','ALPL_gex','VNN2_gex','HLA-DQA1_gex','HLA-DQB1_gex','MARCO_gex','MERTK_gex',
                                     'FCN1_gex','CD2','CD1c','CD1C_gex','S100A9_gex','C1QB_gex','S100A8_gex','C1QA_gex','C1QC_gex','CPA3_gex','C1QTNF4_gex','CD352']) 
h.add_subset('lymphocyte', 'gross', dict(neg=['CD163','MRC1_gex','CD64','OLR1_gex','MPO_gex','ELANE_gex','LAMP3_gex','CLEC9A_gex','CCL19_gex','TPSB2_gex','CPA3_gex','CLEC12A_gex'], 
                                         any_of=['CD3','CD19','CD20','JCHAIN_gex','CD2_gex','LILRA4_gex','PLD4_gex','CD2','CD127','CD352'], n=2))
h.add_subset('myelocyte', 'gross', dict(any_of = [['CD163','CD64','OLR1_gex'],['LAMP3_gex','CLEC9A_gex','CCL19_gex','CD1c'],['HLA-DQA1_gex','HLA-DQB1_gex'],
                                                  ['C1QA','C1QB','C1QC'],['S100A9','S100A8','FCN1'],['MARCO','MRC1','SELENOP','MERTK'],['MPO','ELANE','PRSS57','CYTL1','C1QTNF4']], n=[1,2,2,3,2,3,1], 
                                        neg=['CD2_gex','CD19','CD20','CD127','TPSB2_gex','CD3','MMP9_gex','CPA3_gex','CD2','JCHAIN_gex'], any_ofs_connector='|'))
h.add_subset('mast_cell', 'gross',dict(neg=['MRC1_gex','OLR1_gex','CD64','CLEC9A_gex','LILRA4_gex','PLD4','JCHAIN_gex','CD3','CD2_gex','CD19','CD20','CD127','MPO_gex','ELANE_gex','CD163'],
                                         pos=['CD33'],any_of=['TPSB2_gex','CPA3_gex']))
h.add_subset('neutrophil', 'gross',dict(neg=['MRC1_gex','OLR1_gex','CLEC9A_gex','LILRA4_gex','PLD4_gex','JCHAIN_gex','CD3','CD2_gex','CD19','CD20','CD127','MPO_gex','ELANE_gex'],
                                        any_of=[['S100A8_gex','S100A9_gex'],['MMP9_gex','G0S2_gex','ALPL_gex','VNN2_gex']],n=[1,2])) 

h.add_classification('myeloid','myelocyte',['G0S2_gex','CD1C_gex','S100A9_gex','S100A8_gex','CD14_gex','C1QA_gex','C1QB_gex','C1QC_gex','MERTK_gex','MARCO_gex','MRC1_gex','SELENOP_gex','CLEC9A_gex','CCL19_gex',
                                            'CD64','HLA-DQA1_gex','HLA-DQB1_gex','CD1c','FCN1_gex','MPO_gex','ELANE_gex','PRSS57_gex','CYTL1_gex','LAMP3_gex','C1QTNF4_gex','LILRA4_gex','PLD4_gex',
                                            'JCHAIN_gex','CD123','CLEC12A_gex', 'FCERIA_gex','MS4A7_gex','CD109_gex','BIRC3_gex','STK4_gex','TM4SF1_gex'],clf_kwargs=dict(max_features=.1))
h.add_subset('mo_mac','myeloid',dict(any_of=[['MARCO_gex','MRC1_gex','SELENOP_gex','MERTK_gex'],['C1QA_gex','C1QB_gex','C1QC_gex'],['S100A9_gex','S100A8_gex','CD14_gex','CD64_gex'],['FCN1_gex','MS4A7_gex']],n=[1,2,2,1],
                                     any_ofs_connector='|', neg=['CD1c','CD1C_gex','MPO_gex','ELANE_gex','PRSS57_gex','CYTL1_gex','CCL19_gex','LAMP3_gex','FCERIA_gex','BIRC3_gex','STK4_gex','TM4SF1_gex','CD109_gex']))
h.add_subset('dc','myeloid',dict(any_of=['CD1C_gex','CD1c','CLEC9A_gex','CCL19_gex','LAMP3_gex','FCERIA_gex','CD109_gex','BIRC3_gex','STK4_gex','TM4SF1_gex'],
                                 neg=['MPO_gex','ELANE_gex','PRSS57_gex','CYTL1_gex','FCN1_gex', 'MARCO_gex','SELENOP_gex','MERTK_gex','S100A9_gex','S100A8_gex','MS4A7_gex'],n=1, any_ofs_connector='|'))
h.add_subset('mpdc','myeloid',dict(neg=['MARCO_gex','MRC1_gex','SELENOP_gex','C1QA_gex','C1QB_gex','C1QC_gex','CD14_gex','CD64','FCN1_gex','MPO_gex',
                                        'ELANE_gex','CD109_gex','BIRC3_gex','STK4_gex'],any_of=['LILRA4_gex','PLD4_gex','JCHAIN_gex','CD123'],n=2))

h.add_classification('mono_and_mac','mo_mac',['C1QA_gex','C1QB_gex','C1QC_gex','MARCO_gex','MERTK_gex','SELENOP_gex','FCGR3A_gex','S100A9_gex','S100A8_gex',
                                              'SELL_gex','CD14','CD33_hi','CX3CR1_gex','FCN1_gex','MS4A7_gex','CD16','CD62L','CD14_gex',
                                              'HLA-DQA1_lo','CLEC12A_lo','CD99','APOC1_gex','APOE_gex','ACP5_gex','CD81_gex','LYVE1_gex'],in_danger_noise_checker='noise only')
h.add_subset('macrophage','mono_and_mac',dict(any_of=[['C1QA_gex','C1QB_gex','C1QC_gex','MARCO_gex','SELENOP_gex','APOC1_gex','APOE_gex','LYVE1_gex'],
                                                      ['HLA-DQA1_lo','ACP5_gex','CD81_gex']],neg=['FCN1_gex','CLEC12A_lo','CD99'],n=[1,2],any_ofs_connector='|')) # MS4A7
h.add_subset('monocyte_classical','mono_and_mac',dict(neg=['MERTK_gex','MARCO_gex','SELENOP_gex','MS4A7_gex','CD16','ACP5_gex','CD81_gex','LYVE1_gex'], 
                                                      any_of=[['S100A9_gex','S100A8_gex'],['SELL_gex','CD62L','CLEC12A_lo','CD99','CD14','CD14_gex']],pos=['FCN1_gex','CD33_hi'],n=[2,2],any_ofs_connector='|'))
h.add_subset('monocyte_nonclassical','mono_and_mac',dict(neg=['MARCO_gex','SELENOP_gex','SELL_gex','APOE_gex','APOC1_gex','CD14','LYVE1_gex'],
                                                         any_of=[['CX3CR1_gex','MS4A7_gex','C1QA_gex','C1QB_gex','CD16_gex','FCGR3A_gex'],['FCN1_gex','CLEC12A_lo']],n=[2,1],any_ofs_connector='&'))

h.add_classification("dcs","dc",["XCR1_gex","CLEC9A_gex", "CADM1_gex", "CLEC10A_gex","FCER1A_gex","CD1C_gex","CCR7","CCL19_gex"])
h.add_subset("dc1","dcs",dict(any_of=["XCR1_gex","CLEC9A_gex","CADM1_gex",], neg=["CLEC10A","FCER1A","CD1C_gex"]))
h.add_subset("dc2","dcs",dict(any_of=["CLEC10A_gex","FCER1A_gex","CD1C_gex"], neg=["XCR1_gex","CLEC9A_gex","CADM1_gex"]))
h.add_subset("dc_migratory","dcs",dict(any_of=["CCR7","CCL19_gex"], neg=["XCR1_gex","CLEC9A_gex","CADM1_gex"]))

h.add_classification('lymphoid','lymphocyte',['CD3','TCR_ab','CD19','CD20','MZB1_gex','JCHAIN_gex','CD2_gex','TRDV1_gex','TRDV2_gex','TRDV3_gex','KLRF1_gex','IL7R_gex','KIT_gex','IL4I1_gex','RORC_gex',
                                              'NCR2_gex','LILRA4_gex','PLD4_gex','TCR_Vd2'])
h.add_subset('t_cell', 'lymphoid', dict(neg = ['CD19','CD20','MZB1_gex','JCHAIN_gex','NCR2_gex','LILRA4_gex','PLD4_gex'], 
                                        any_of = [['CD3'],['TRDV1_gex', 'TRDV2_gex', 'TRDV3_gex']],any_ofs_connector='|')) 
h.add_subset('b_cell_like', 'lymphoid', dict(any_of = ['CD19', 'CD20', 'MZB1_gex', 'JCHAIN_gex','LILRA4_gex','PLD4_gex'], neg=['CD3', 'IL7R_gex'])) 
h.add_subset('nk_ilc', 'lymphoid', dict(neg = ['CD3', 'CD19', 'CD20', 'MZB1_gex', 'TCR_Vd2','TRDV1_gex', 'TRDV2_gex', 'TRDV3_gex', 'JCHAIN_gex','LILRA4_gex','PLD4_gex'], 
                                        any_of=['IL7R_gex', 'KIT_gex', 'IL4I1_gex', 'RORC_gex', 'KLRF1_gex', 'NCR2_gex']))

h.add_classification('tcr','t_cell',['TCR_ab','TRAC_gex','TRDC_gex','TCR_Vd2','TRDV1_gex', 'TRDV2_gex', 'TRDV3_gex'],in_danger_noise_checker=False, clf_kwargs=dict(max_features=.25))
h.add_subset('ab_t','tcr',dict(any_of=['TCR_ab','TRAC_gex'], neg=['TRDC_gex','TCR_Vd2','TRDV1_gex', 'TRDV2_gex', 'TRDV3_gex']))
h.add_subset('gd_t','tcr',dict(neg=['TCR_ab'], any_of=['TRDC_gex','TCR_Vd2','TRDV1_gex', 'TRDV2_gex', 'TRDV3_gex'],n=2))

h.add_classification('cd4cd8','ab_t',['CD4','CD8','CD4_gex','CD8A_gex'], in_danger_noise_checker=False, clf_kwargs=dict(max_features=.1)) 
h.add_subset('cd4_t','cd4cd8',dict(any_of=['CD4','CD4_gex'],neg=['CD8','CD8A_gex']))
h.add_subset('cd8_t','cd4cd8',dict(any_of=['CD8','CD8A_gex'],neg=['CD4','CD4_gex']))

h.add_classification('cd4_mem','cd4_t',['FOXP3_gex','CD62L','CD45RA','CCL5_gex','SELL_gex','CD127','CCR7_gex','CD69',"CTLA4_gex","GZMB_gex","GZMK_gex","RORC_gex","FCGR3A_gex","GZMA_gex","IFNG_gex","IKZF2_gex",
                                        "IL10_gex","PDCD1_gex","FGFBP2_gex","GZMH_gex","CCR7","KLRG1_gex"], clf_kwargs=dict(max_features=.1))
h.add_subset('cd4_naive_cm','cd4_mem',dict(any_of=['CD62L',"SELL_gex","CCR7"],neg=['FOXP3_gex', "IKZF2_gex", "CCL5_gex", "FCGR3A_gex", "FGFBP2_gex", "KLRG1_gex",])) 
h.add_subset('cd4_treg','cd4_mem',dict(any_of=[["FOXP3_gex", "IL10_gex", "IKZF2_gex",]], neg=["CCL5_gex", "CD127"], n=[1]))
h.add_subset('cd4_em','cd4_mem',dict(any_of=['CD69', "GZMA_gex", "GZMB_gex", "GZMK_gex", "GZMH_gex", "IFNG_gex", "PDCD1_gex"], neg=['CD62L', 'CD45RA', 'FOXP3_gex', "IKZF2_gex"]))
h.add_subset('cd4_temra','cd4_mem',dict( any_of=[['CD45RA',],['CCL5_gex', "FCGR3A_gex", "FGFBP2_gex", "KLRG1_gex", "GZMA_gex", "GZMB_gex", "GZMK_gex", "GZMH_gex", "IFNG_gex", "PDCD1_gex"]],
                                         neg=['FOXP3_gex', 'CD62L', 'SELL_gex', 'CCR7_gex', "IKZF2_gex"], n=[1,1], any_ofs_connector='&'))

h.add_classification('cd4_ncm', 'cd4_naive_cm', ['CD45RA',"CD95"])
h.add_subset('cd4_naive','cd4_ncm',dict(pos=['CD45RA'], neg=["CD95"]))
h.add_subset('cd4_cm','cd4_ncm',dict(neg=['CD45RA'],pos=["CD95"]))

h.add_classification('cd4_res','cd4_em',['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6_gex'])
h.add_subset('cd4_trm','cd4_res',dict(any_of=['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6_gex'],n=2))
h.add_subset('cd4_tem','cd4_res',dict(neg=['CD69','CD103','CD49a','ITGA1_gex','CXCR6_gex']))

h.add_classification('cd8_mem','cd8_t', ['CD62L','CD45RA',"KLRG1_gex","CD103",'CD57','CCL5_gex','SELL_gex','CCR7_gex',"CD69",
                                         "GZMB_gex","GZMK_gex","GZMA_gex","IFNG_gex","PDCD1_gex","FCGR3A_gex","NKG7_gex","FGFBP2_gex"],
                                         in_danger_noise_checker=False, clf_kwargs=dict(max_features=.1)) 
h.add_subset('cd8_naive_cm','cd8_mem',dict(pos=['CD62L'], neg=['CCL5_gex',"NKG7_gex", "FCGR3A_gex",'CD57',"FGFBP2_gex","KLRG1_gex"]))
h.add_subset('cd8_tem_mait','cd8_mem',dict(any_of=['CD69','CCL5_gex',"GZMB_gex","GZMK_gex","GZMA_gex","IFNG_gex","PDCD1_gex","KLRG1_gex","CD57","CD103"], neg=['CD62L','CD45RA'],n=1))
h.add_subset('cd8_temra','cd8_mem',dict(any_of=[['CD45RA'],['CCL5_gex',"FCGR3A_gex","GZMB_gex","GZMK_gex","GZMA_gex","IFNG_gex","PDCD1_gex",'CD57',"FGFBP2_gex","KLRG1_gex"]],
                                        neg=['CD62L','SELL_gex',"CCR7_gex"], n=[1,1], any_ofs_connector='&'))

h.add_classification('cd8_ncm','cd8_naive_cm', ['CD45RA', "CD95", "CD45RO", "CCL5_gex"])
h.add_subset('cd8_naive','cd8_ncm', dict(pos=['CD45RA'], neg=["CD95","CD45RO","CCL5_gex"]))
h.add_subset('cd8_cm','cd8_ncm', dict(any_of=["CD95","CD45RO","CCL5_gex"], neg=['CD45RA']))

h.add_classification('mait','cd8_tem_mait',['TRAV1-2_gex','DPP4_gex','SLC4A10_gex','KLRB1_gex','TCR_Va7_2','TRAC_gex','TCR_ab'], clf_kwargs=dict(max_features=.1))
h.add_subset('cd8_mait','mait',dict(any_of=['TRAV1-2_gex','DPP4_gex','SLC4A10_gex','KLRB1_gex','TCR_Va7_2'], n=3))
h.add_subset('cd8_not_mait','mait',dict(neg=['TRAV1-2_gex','DPP4_gex','SLC4A10_gex','KLRB1_gex','TCR_Va7_2'], any_of=['TRAC_gex','TCR_ab'], n=1))

h.add_classification('cd8_res','cd8_not_mait',['CD69','ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6_gex'])
h.add_subset('cd8_trm','cd8_res',dict(any_of=[['CD69'],['ITGAE_gex','CD103','ITGA1_gex','CD49a','CXCR6_gex']],n=[1,2],any_ofs_connector='|'))
h.add_subset('cd8_tem','cd8_res',dict(neg=['CD69','CD103','CD49a','ITGA1_gex','CXCR6_gex']))

h.add_classification('b_c_like','b_cell_like',['CD19','CD20','JCHAIN_gex','JCHAIN_hi','MZB1_gex','LILRA4_gex','PLD4_gex','CD319_CRACC','CD123','CLEC12A_gex'])
h.add_subset('b_cell','b_c_like',dict(neg=['MZB1_gex','JCHAIN_gex','LILRA4_gex','CD319_CRACC','CD123','CLEC12A_gex'],any_of=['CD19','CD20']))
h.add_subset('plasma','b_c_like',dict(any_of=['JCHAIN_hi','MZB1_gex','JCHAIN_gex','CD319_CRACC'],neg=['LILRA4_gex','PLD4_gex','CD20','CD123','CLEC12A_gex'],n=2))
h.add_subset('pdc','b_c_like',dict(neg=['CD19','CD20','JCHAIN_hi'],any_of=['LILRA4_gex','PLD4_gex','JCHAIN_gex','CD123','CLEC12A_gex'],n=2))

h.add_classification('b_mem','b_cell', ['IgD','CD27','AICDA_gex',"SUGCT_gex","RGS13_gex","LMO2_gex","ELL3_gex","CD38",'TCL1A_gex','YBX3_gex',
                                        'IL4R_gex','ITGAX_gex','TBX21_gex','CD11c', 'FGR_gex','IGHG1_gex','IGHA1_gex','IGHG2_gex','IGHA2_gex',
                                        'CD27_gex','SCIMP_gex','BACH2_gex','PLPP5_gex'], in_danger_noise_checker='noise only')
h.add_subset('b_naive','b_mem',dict(pos=['IgD'], any_of = ['TCL1A_gex','YBX3_gex','IL4R_gex','BACH2_gex','PLPP5_gex'], neg=['CD27','CD27_gex','SCIMP_gex','IGHG1_gex','IGHA1_gex','IGHG2_gex','IGHA2_gex'],n=1))
h.add_subset('b_memory','b_mem',dict(pos=['CD27_gex'],any_of=['CD27','SCIMP_gex','IGHG1_gex','IGHA1_gex','IGHG2_gex','IGHA2_gex'],n=2))
h.add_subset('b_gc','b_mem',dict(any_of=["AICDA_gex","SUGCT_gex","RGS13_gex","LMO2_gex","ELL3_gex",],n=2))
h.add_subset('b_age','b_mem',dict(any_of=["CD11c",'ITGAX_gex','TBX21_gex','FGR_gex'],n=2))

h.add_classification("plasma_blast","plasma",["SDC1_gex","MKI67_gex","TOP2A_gex","TNFRSF18_gex"])
h.add_subset("plasma_cell","plasma_blast", dict(any_of=["SDC1_gex","TNFRSF18_gex"],neg=['MKI67_gex',"TOP2A_gex"]))
h.add_subset("plasmablast","plasma_blast", dict(any_of=['MKI67_gex',"TOP2A_gex"],neg=["SDC1_gex"]))

h.add_classification('nk_ilc_subsets','nk_ilc',['EOMES_gex','GZMH_gex', 'IL7R_gex', 'FCGR3A_gex', 'GZMK_gex', 'GZMB_gex', 'KIT_gex','IL4I1_gex','RORC_gex','NCR2_gex','CD335_NKp46','CD103',"CD57","SELL_gex","CD62L","CD69"])
h.add_subset('nk_cd56dim','nk_ilc_subsets',dict(pos=['FCGR3A_gex','GZMB_gex','EOMES_gex'], neg=['GZMK_gex','KIT_gex','IL4I1_gex','RORC_gex','NCR2_gex']))
h.add_subset('nk_cd56br','nk_ilc_subsets',dict(any_of=['GZMK_gex','EOMES_gex'], pos=['CD335_NKp46'], neg=['FCGR3A_gex','KIT_gex','IL4I1_gex','RORC_gex','NCR2_gex','CD103']))
h.add_subset('ilc1','nk_ilc_subsets',dict(any_of=['NCR2_gex'], neg=['FCGR3A_gex' ,'EOMES_gex','KIT_gex']))  
h.add_subset('ilc3','nk_ilc_subsets', dict(any_of=['IL7R_gex','KIT_gex','IL4I1_gex','RORC_gex'], neg=['FCGR3A_gex', 'GZMK_gex', 'GZMH_gex', 'GZMB_gex','EOMES_gex','NCR2_gex'],n=2))
h.add_subset('nk_ilc_precursor','nk_ilc_subsets', dict(any_of=['SELL_gex','CD62L'],neg=['FCGR3A_gex','NCR2_gex',"RORC_gex","CD57","CD69"]))
```

## Prior Versions

[Placeholder]