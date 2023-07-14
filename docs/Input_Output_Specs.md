# Input and Output Specifications
## Expected Inputs

<img align="left" src="https://raw.githubusercontent.com/scverse/anndata/main/docs/_static/img/anndata_schema.svg" width="250" style="border:4px solid black;background-color:white">

### General input structure
MMoCHi expects an `anndata.AnnData` object ([AnnData](https://anndata.readthedocs.io/en/latest/)), and requires ~3-5 times the amount of RAM needed to store this object in memory to run. Currently for multimodal CITE-Seq classification, MMoCHi expects the `.X` of the `AnnData` object to contain gene expression (GEX) data and the `.obsm[data_key]` to contain a `pandas.DataFrame` ([pandas](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)) with protein expression (ADT) data. This expression data should be library-size normalized and batch corrected, as necessary.

Helper functions (such as `mmc.preprocess_adatas`) may help convert `AnnData` objects into the right format, although you should check to make sure the object is correctly formatted and data are correctly transformed. 

### `.X`
- In the `.X` of the `AnnData`, MMoCHi classification requires library-size normalized feature expression
    - Most often this would be log-normalized gene expression stored in a `scipy.sparse_matrix` but could also be stored in a `numpy.array`

### `.obsm`
- Expression of the second modality can be stored in the `.obsm[data_key]` as a `pandas.DataFrame` ([pandas](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)).
    - Most often this would be log-normalized protein expression. Protein expression can also be batch-corrected, for example by using `mmc.landmark_register()`
    - The DataFrame should have column names corresponding to features and be indexed identically to `anndata.AnnData.obs_names`. Both must be unique and encoded as strings.
    - For convenience, the default `data_key` for most functions can be edited by setting `mmc.DATA_KEY='new_data_key'`

### `.obs`
- The index of the `.obs` (also stored in `anndata.AnnData.obs_names`) should correspond to cell-barcode identity. These must be unique and encoded as strings.
    - For `.obs_names` that are not unique, this can be fixed using `adata.obs_names_make_unique()`
    - For `.obs_names` that are encoded as integers, this can be changed using `adata.obs_names = adata.obs_names.astype(str)`
    
- For integration of protein expression by landmark registration, or batch-integrated classification, MMoCHi expects a column in the `.obs` delineating batch.
    - This column name is specified by the `batch_key` argument across all MMoCHi functions. Potentially useful batch keys may be donor_id, sequencing_technology, or sample_type.
    - For convenience, the default `batch_key` for most functions can be edited by setting `mmc.BATCH_KEY='new_batch_key'`
    
- During high-confidence thresholding and classification, MMoCHi has built in opportunities to compare classifications to a reference column in the `.obs`. 
    - This is entirely optional but can be useful for troubleshooting. Often good reference columns are cursory manual annotations, clustering, or sample metadata.
    
- In the hierarchy the user can specify cell metadata (such as tissue-of-origin) that can be used to select events for high-confidence populations. 
    - This metadata should be included in the `.obs` and be encoded as strings. An example can be found in the [Hierarchy Design tutorial](Hierarchy_Design.ipynb).
    
### `.var`
- The index of the `.var` (also stored in `anndata.AnnData.var_names`) should correspond to feature_names. These must be unique and encoded as strings.
    - For `.var_names` that are not unique, this can be fixed using `adata.var_names_make_unique()`
    - For `.var_names` that are encoded as integers, this can be changed using `adata.var_names = adata.var_names.astype(str)`

- If you wish to filter out specific features from the `.X` before training the classifier, a boolean mask can be included as a column in the `.var`.
    - Consider filtering out non-protein coding genes (if their expression is noisy and not useful for cell-type discrimination) or known sequencing artifacts.
    
- There is experimental support for expression of multiple modalities to be stored in the `.X`. This may not run reliably. 
    - For this, a column in the `.var`, specifying modality (e.g. `features_type`), must be included.
 
### `mmc.Hierarchy`
 - A `mmc.Hierarchy` object must be created as detailed in the docs, [Integrated Classification tutorial](Integrated_Classification.ipynb), and [Hierarchy Design tutorial](Hierarchy_Design.ipynb).
     - The structure and high-confidence markers used in the hierarchy can be checked using its `.display()` method.
 
## Outputs

### `.obsm['lin']` (after running `mmc.classify()`)

- Named `lin` by default (for lineage), this is a `pandas.DataFrame`, indexed by `anndata.AnnData.obs_names`. 
    - At each level of the hierarchy, columns will be added to this DataFrame containing important results and information.
    - After training a classifier, the `.obsm['lin']` will include columns for: 
        - `{level}_hc`: String, High-confidence threshold identity (either a cell type, '?' (for events not identified as a high-confidence cell type), or 'nan' for cells not part of the parent subset)
        - `{level}_holdout`: Bool, True for events that are *explicitly* held out from training (i.e. not used for random forest training, nor training data selection)
        - `{level}_train`: Bool, True for events that are used for random forest training
        - `{level}_traincounts`: Integer, The number of times an event was duplicated within the training dataset.
    - After using the classifier for prediction, the `.obsm['lin']` will include columns for: 
        - `{level}_class`: String, The classification at the given level, or 'nan' for cells not part of the parent subset.
        - `{level}_probability`: Float, The classifier's confidence of the predicted subset being correct for a given event.
    - Note, cutoff nodes (e.g. when using `h.add_classification(is_cutoff=True)`) in the hierarchy will only produce `{level}_hc` and `{level}_class` columns when run.
    - Events that are identified as 'noise' or excluded during subsampling during training data selection will be `False` for both `{level}_holdout` and `{level}_train`

### `.obs` (after running `mmc.terminal_names()`)

- The `classification` column is the lowest-level annotation of a cell type. If two annotations are put in parallel on the hierarchy, it will correspond to the annotation specified last.
    
- The `certainty` column corresponds to the Random Forest's confidence level (scaled 0-1) for the classified subset identity at the most specific level of annotation. 
    - Without calibration enabled, this corresponds to the percent of the forest in agreement. With calibration, this should more-closely match the certainty of the given cell call.
    - This will be higher for cells that are well-represented by high-confidence events and lower for cells that are less well-represented and may be useful for identifying problem-areas in classification
    - This should not be as a direct readout of classification certainty, as it does not take into account certainty of intermediate calls in the hierarchy. The `probability_cutoff` parameter in `mmc.classify()` may provide more ability to select for cells classified with high certainty. 

### `.uns[{level}_probabilities]` (after running `mmc.classify(probabilities_suffix='_probabilites')`)

- MMoCHi will optionally add `pandas.DataFrame` objects encoding the predicted probabilities for all classes in the multiclass at each level to the `.uns`. 
    - This can be very useful for troubleshooting classification results or using MMoCHi to generate cell type scores.

### `mmc.Hierarchy` (after running `mmc.classify()`)

 - As before, the `mmc.Hierarchy` object will contain information about each classification level and its child nodes.
 - While running `mmc.classify()`, the features used, random forest classifier, and any calibration will be saved into the hierarchy object. These can be accessed to identify important features for classification, as demonstrated in the [Feature importance tutorial](Exploring_Feature_Importances.ipynb).
 - The hierarchy object produced by `mmc.classify()` has all items needed to apply this classification on held-out data or other datasets, by running `mmc.classify(retrain=False)`, as in the [Pretrained Classification tutorial](Pretrained_Classification.ipynb). Projecting the classification requires that expression information for all features used in the original classification are available in the held-out dataset.
    - One should also be cautious to ensure that cell types present in the training dataset are representative of cell types present in the held-out data. 
    - Rather than applying pretrained MMoCHi models, we recommend attempting to retrain new classifiers at each level of the hierarchy using high-confidence thresholding on this new dataset. 
 
### `.log` file (if `mmc.log_to_file()` is run before `mmc.classify()`)

-  While functions in the `mmochi` package are running, they will use logging to print info, warnings, and errors to the output or error stream. If you run `mmc.log_to_file()` these outputs will also be written to that file, along with other information that is not printed to the output stream. 

## Example code for formatting 

To help with formatting your data properly for MMoCHi (with one modality stored in the `.obsm`) we have provided some example cases below.

If you have two separate `anndata.AnnData` objects indexed identically by cell barcode, you can combine them as follows:
```python
import scanpy as sc
import pandas as pd

adata = sc.read_adata('gex.h5ad')
protein=sc.read_adata('adt.h5ad')
# Note if protein.X is a sparse matrix, you may need to use protein.X.A
adata.obsm['protein'] = pd.DataFrame(protein.X, index=protein.obs_names, columns = protein.var_names)
adata.write_h5ad('mmochi_formatted.h5ad')
```

If you have two separate `.loom` files, indexed identically by cell barcode, you can follow the following code, replacing `ra.Antigen` with the row attribute containing protein names:
```python
import scanpy as sc
import loompy 
import numpy as np
import pandas as pd

adata = sc.read_loom('gex.loom')
with loompy.connect('adt.loom') as protein:
    # This can be helpful to remove '/' characters from antigen names.
    proteins = [p.replace('/','_') for p in protein.ra.Antigen]
    adata.obsm['protein'] = pd.DataFrame(protein[:,:], 
                                         index=adata.obs_names, 
                                         columns=proteins)
adata.write_h5ad('mmochi_formatted.h5ad')
```

If you have a `mudata.MuData` object, you can follow the following code, replacing the `'gex'` and `'adt'` modalities with your modality names:
```python
import mudata as md
mdata = md.read_h5mu('mudata.h5mu')

adata = mdata.mod['gex'].copy()
# Note if mdata.mod['adt'].X is a sparse matrix, you may need to use protein.X.A
adata.obsm['protein'] = pd.DataFrame(mdata.mod['adt'].X, 
                                     index=mdata.mod['adt'].obs_names, 
                                     columns=mdata.mod['adt'].var_names)
adata.write_h5ad('mmochi_formatted.h5ad')
```