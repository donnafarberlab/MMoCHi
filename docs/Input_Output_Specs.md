# MMoCHi
# Expected Inputs

<img align="left" src="https://raw.githubusercontent.com/scverse/anndata/main/docs/_static/img/anndata_schema.svg" width="250" style="border:4px solid black;background-color:white">

## General input structure
MMoCHi expects an `anndata.AnnData` object ([anndata](https://anndata.readthedocs.io/en/latest/)), loaded into memory. During the course of running MMoCHi training or prediction, there will need to be sufficient memory available to duplicate all features used.

All inputted expression to MMoCHi should be library-size normalized, and batch corrected, as necessary.

Currently for multimodal CITE-Seq classifiation, MMoCHi expects the `.X` of the `AnnData` object to contain gene expression, and a `pandas.DataFrame` ([pandas](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)) in the `.obsm[data_key]` that contains protein expression.

Much of the formatting required can be automated using helper functions (such as `mmc.preprocess_adatas`), although you should manually check to make sure the object is correctly formatted and all data are correctly transformed. 

# Detailed specification of inputs

#### `.X`
- In the `.X` of the `anndata`, MMoCHi classification requires library-size normalized feature expression
    - Most often this would be log-normalized gene expression stored in a `scipy.sparse_matrix`, but could also be stored in a `numpy.array`

#### `.obsm`
- Expression of the other modalities can be stored in  `pandas.DataFrame` ([pandas](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)) objects in the `.obsm[data_key]`. MMoCHi currently only supports 1 modality stored in the `.obsm`
    - Most often this would be log-normalized protein expression. Protein expression can also be batch-corrected using `mmc.landmark_register()`
    - Column names must correspond to the feature names, and row names should correspond to the `anndata.AnnData.obs_names`

#### `.obs`
- The index of the `.obs`, also stored in `anndata.AnnData.obs_names`, should correspond to cell-barcode identity. While the actual values of these are not important, they must be unique and a string, not integers.
    - For `.obs_names` that are already encoded as integers, this can be changed using `adata.obs_names = adata.obs_names.astype(str)`
    
- For landmark registration of protein expression or batch-integrated classification, MMoCHi expects a column in the `.obs` delineating batch.
    - This column name is specified by the `batch_key` argument across all MMoCHi functions. Potentially useful batch keys may be donor_id, sequencing_technology, or sample_type.
    
- During ground truth-thresholding and classification, MMoCHi has built in steps to make comparisons between the classifications and a reference column in the `.obs`. 
    - This is entirely optional, but can be very useful for troubleshooting classifications. Often good reference columns are cursory manual annotations, clustering, or sample metadata.
    
- In the hierarchy the user can specify cell metadata (such as tissue-of-origin) that can be used to select events for ground-truth. These columns work best when encoded as strings, and must be included in the `.obs`.
    
#### `.var`
- The index of the `.var`, also stored in `anndata.AnnData.var_names`, should correspond to feature_names. While the actual values of these are not important, they must be unique and a string, not integers.

- If specific features should be filtered out of the `.X` for training the classifier, a mask column of features to include can be provided. 
    - Some genes to consider filtering out of MMoCHi training include non-protein coding genes (if their expression is noisy and not useful for cell-type discrimination) or known sequencing artifacts).
    
- Expression of multiple modalities can be stored in the `.X`. When this is the case, a column in the `.var`, specifying modality (e.g. `features_type`), must be included.
    - Note this feature is currently experimental, and may not function through every step of classification.
 
#### `mmc.Hierarchy`
 - The `mmc.Hierarchy` object must be created, as specified in the docs and [tutorials](/docs/Classifier_Demo.ipynb).
     - The structure and ground-truth markers used in the hierarchy can be checked using its `.display()` method.
 
# Specification of Outputs

#### `.obsm['lin']` (after running `mmc.classify()`)

- Named `lin` by default (for lineage), this is a `pandas.DataFrame`, indexed by `anndata.AnnData.obs_names`, 

#### `.obs` (after running `mmc.terminal_names()`)

- The `classification` column is the lowest-level annotation of a celltype. If two annotations are put in parallel on the hierarchy, it will correspond to the annotation specified last.
    
- The `conf` column corresponds to the Random Forest's confidence level (scaled 0-1) at the lowest-level of annotation. 
    - This will be higher for cells that are well-represented by ground-truth and lower for cells that are less well-represented, and may be useful for identifying problem-areas in classification
    - This should not be used for direct annotation of cells as "high-confidence" as it does not take into account confidence scores elsewhere on the hierarchy. More success for "high-confidence" annotations may come from using the `probability_cutoff` parameter in `mmc.classify()`

#### `.uns[{...}_proba]` (after running `mmc.classify()` with `proba_suffix` not equal to `None`)

- MMoCHi will optionally add `pandas.DataFrame` objects to the unstructured for each classification level, corresponding to the predicted probabilities for all classes in the multiclass at each level. 
    - This can be very useful for troubleshooting classification results, or using MMoCHi to "score" cell types based on an identity.

#### `mmc.Hierarchy` (after running `mmc.classify()`)

 - As before, the `mmc.Hierarchy` object will contain information about each classification level and its child nodes.
 - While running `mmc.classify()`, the features used, random forest classifier, and any calibration will be saved into the hierarchy object. These can be accessed to identify important features for classification, as demonstrated in [tutorials](/docs/Feature_Importances.ipynb).
    - The hierarchy object produced by `mmc.classify()` has all items needed to apply this classification on held-out data, or other datasets, by using the parameter `retrain = False` on `mmc.classify()`. Projecting the classification requires that expression information for all features used in the original classification are available in the held-out dataset.
    - One should also be cautious to ensure that celltypes present in the training dataset are representative of celltypes present in the held-out data. If applying to equivalent CITE-Seq data, one can also rerun `mmc.classify()` with `retrain = True` to retrain new classifiers at each level of the hierarchy. 
 
#### `.log` file (if `mmc.log_to_file()` is run before `mmc.classify()`)

-  While functions in the `mmochi` package are running, they will use logging to print info, warnings, and errors to the output or error stream. If you run `mmc.log_to_file()` these outputs will also be written to that file, along with other information that is not printed to the output stream. 