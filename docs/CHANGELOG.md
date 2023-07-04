# Change Log

All notable changes to this package will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/) and after v1.0.0 release, 
this project will adhere to [Semantic Versioning](http://semver.org/).

To update, replace your repository with a fresh clone, and install using pip, as before:

```
git clone https://github.com/donnafarberlab/mmochi.git
conda activate mmochi
cd mmochi
pip install .
```
---
## Current version

### [0.2.2] - 30JUN23

#### Added

- Automated documentation hosted by ReadTheDocs

- Issue templates

- Examples MMoCHi hierarchies

- Examples in Input/Output Specifications for how to correctly format various objects

#### Fixed 

- Updated documentation of various functions

- Added temporary requirement limiting scikit-learn to below 1.3.0, as that update breaks imbalanced-learn

- Removed mmochi package from Python3_8_requirements.txt 

- Prevented plot windows from oppening while running pytest on some systems

### [0.2.1] - 17JUN23

#### Added

- Functionality to restore `adata.obs_names` if `mmc.classify` is interrupted. 
    - Internally, `mmc.classify` converts `adata.obs_names` to indicies, and remaps them at the end of the function.
    - Now, `adata.obs_names` can be restored from the temporary column: `adata.obs[‘MMoCHi_obs_names’]`
    
- New option for `mmc.density_plot` whether to hide events with 0 expression (default: True)

#### Changed

- Nomenclature changed for columns in the `.obsm`: the `'_tcounts'`  to `'_traincounts'`, `'_proba'` to `'_probability'` or `'_probabilities'` (in `.uns`), `'conf'` default column name to `'certainty'` (in the `.obs`). Note that although MMoCHi functions will not support this old format, previously generated classifications can be renamed to match this new format.

- Default log-normalization in preprocess_adatas changed for GEX to counts per 10k and ADT to counts per 1k

- Updated tutorials, doctstrings, README, and input_output_specs for clarity

#### Fixed

- Removed `leidenalg` from `requirements.txt` to prevent conda install from hanging while trying to solve the environment

- Bug handling `np.nan` in `mmc.utils.umap_interrogate_level`

- Disabled some internal warnings in `mmc.utils`

### [0.2.0] - 26MAY23

#### Added

- Introduced new tutorials for [Hierarchy Design](Hierarchy_Design.ipynb), [High Confidence Thresholding](High_Confidence_Thresholding.ipynb), [Exploring Feature Importances](Exploring_Feature_Importances.ipynb), and [Pretrained Classification](Pretrained_Classification.ipynb)

- Edited the tutorials for [Integrated Classification](Integrated_Classification.ipynb) and [Landmark Registration](Landmark_Registration.ipynb)

- Created a new helper function (`mmc.umap_interrograte_level()`) for plotting UMAPs to understand classification performance

- Added variables for setting the default `data_key` and `batch_key` for many functions by changing `mmc.DATA_KEY` and `mmc.BATCH_KEY` for many functions. Note, defaults for these two for some functions have now been changed

- Added checks to ensure `_obs`, `_gex`, and `_mod_` are not in the feature names inputted into MMoCHi

#### Changed

- 'Ground Truth' and 'gt' have been renamed to 'High Confidence' and 'hc' throughout the package for clarity

- 'hold_out' was replaced with 'holdout' for consistency

- Nomenclature of column names in `adata.obsm['lin']` has changed. Note that although MMoCHi functions will not support this old format, previously generated .obsm['lin'] dataframes contain enough data to be converted to the new format.
    - Columns in the `.obsm['lin']` now include (see [Input/Output Specifications](Input_Output_Specs.md) for more detail): 
        - `{level}_hc`: High-confidence threshold identity
        - `{level}_holdout`: Events *explicitly* held out from training data selection
        - `{level}_train`: Events used for random forest training
        - `{level}_tcounts`: Number of times events were duplicated during random forest training.
    - Notably, events that are identified as 'noise' or excluded during subsampling during training data selection will be `False` for both `{level}_holdout` and `{level}_train`

- Coloring of `mmc.run_all_thresholds()` histograms was reversed to match those of `mmc.umap_thresh()`.

- Made many functions private (this should not affect any functions the user interacts with).

#### Fixed

- Events held out from training were sometimes given `NaN` for the `{level}_tcounts` and `{level}_train`, which could be mistaken when typecasting to bool

- Added garbage collection to reduce wasteful memory usage buildup after plotting UMAPs (due to occasional object duplication).

- Improved pytesting to validate post-classification structure of the `.obsm['lin']`

- Fixed typos in some log and print statements

- Updated docstrings and tutorials for compatability with sphinx.

- Replaced explicit reference to `.obsm['lin']` with `key_added` parameter for `mmc.plot_important_features()`

#### Removed

- Unneccesary call to `clf.feature_importances_` in `mmc.feature_importances()`

- Removed `load_dir` argument in `mmc.Hierarchy()`. Hierarchies can now be loaded from a full file path provided by the `load` parameter. 

- Removed `optional-requirements.txt`. All mandatory requirements are now listed in `requirements.txt`


### [0.1.4] - 09MAR23

#### Added
 
- Support for predicting on extra-large datasets stored in int64-indexed sparse matrices. (These are not supported by scikit-learn, so they are split into 
  bite-sized chunks.)
  
- Loading saved hierarchies from other directories, using the new `load_dir` argument to the `mmc.Hierarchy()` constructor
  
#### Changed
  
- During ground-truth cleanup, the PCA is now run on ***scaled*** highly variable features, such that highly expressed features, or differences in expression levels
  between modalities do not dominate.
  
- Disable `reduce_features_min_cells` in `mmc.classify` when `retrain == True`, so that features are not filtered out when projecting a classifier onto a new dataset.
  If highly expressed features need to still be removed, this can be performed prior to inputting into the `mmc.classify` function
 
#### Fixed

- `**kwargs` can now be passed through `mmc.plot_confusion` to `sklearn.metrics.ConfusionMatrixDisplay.from_predictions()`

- Added support for recent versions of scikit-fda, which should include support for the current version: `skfda==0.8.1`. 

---

## Past versions

### [0.1.3] - 22Dec22

#### Added
 
- Updated tutorials for [Integrated Classification](Integrated_Classification.ipynb) and [Landmark Registration](Landmark_Registration.ipynb)
  
- Docstrings and typing hints have been written for most user functions, and many internal functions for code clarity.

- The `n_estimators` used for each batch during batch-integration can now be weighted by representation in the dataset (`weight_integration=True`), or kept 
  weighted equally in the forest (old performance, and the default)

- `min_cluster_size` parameter to `mmc.classifier.borderline_balance`

- added `base` parameter to `mmc.hierarchy.Hierarchy.get_clf` function to (if needed) strip the outer calibration classifier, and return a random forest

- Added `mmc.landmarking.update_peak_overrides` as a convenience function for creating a dictionary of manual peak overrides
  , `mmc.landmarking.update_landmark_register` as a function to test new landmark registration settings on a single batch-marker (after landmark registration 
  has been run on the entire dataset), and 
  
- Peak overrides can now be specified as single-positives, by passing a `[None, float]`

- A new option for density plotting, `mmc.landmarking.density_plot_total`, has been added, to display the density of a single batch, single marker, in front of
  the density plot of the rest of the dataset (useful when integrating in one new batch)
  
- Added `mmc.landmarking.save_peak_overrides` and `mmc.landmarking.load_peak_overrides` to save and load this object as a JSON. There is also added support for
  defining peak_overrides in `mmc.landmarking.landmark_register_adts` using the path to a JSON for automatic loading.

#### Changed
  
- Moved most classification settings from `mmc.classifier.classify` to hierarchy to be edited on a per-node basis.

- During classifier batch-integration, the `n_estimators` defined in the `clf_kwargs` now refers to the total number of trees in the forest (for more consistent
  performance with classification without batch-integration).
 
- Removed the option for sigmoid calibration to address inconsistencies with calibration performance and imbalanced classes. Advanced calibration settings will
  be revisited in a future version. Currently, the calibrated classifier is  trained on data that was not used for initial training, was defined in ground 
  truthing, and (assuming there are enough events) only uses 1/3rd of that data.
  
- The keywords for `batch_id` and `batch` have been updated to `batch` and `batch_key`  in `mmc.landmarking.stacked_density_plots` 
  and `mmc.landmarking.density_plot` for consistency with other functions.
  
- Changed the order of returns for `mmc.thresholding.threshold` so that regardless of whether `run=False` or `run=True`, thresh is the first returned object, 
  and the thresholded data is an optional extra returned object.
 
#### Fixed

- Improved pytest coverage of all modules, to test broader use-cases.

- Specified peak overrides in `mmc.landmarking._landmark_registration` should now more accurately reflect the closest possible mapping of those values to the 
  proper position in the FDA function.

- Using `return_graph=True` and `plot=False` together on `mmc.hierarchy.Hierarchy.display` will now return None (expected behavior) instead of throwing an
  error. 

- When initializing a new `mmc.hierarchy.Classification` with a predefined classifier, there is now error checking to ensure feature_names are also provided

- Warnings about any cython import failures during the import of scikit-fda are now silenced.

- Fixed performance of `mmc.utils.generate_exlcusive_features` if `adatas` is given as a list of str. Previously this would mistakenly return an empty list.

- Fixed reading in 10x formatted .h5 files without url backend using `mmc.utils.preprocess_adatas`

- Fixed `get_data` function when passed a keyword including `_obs` in the variable name

- Added fixes for cases where thresholds in `umap_thresh` were out of the bounds of the data.

- Fixed error where if features were limited on a per-classification-level basis, the wrong set of features were passed to generate the training matricies.

- Fixed error if no clusters needed to be balanced in `mmc.classifier.balance_borderline`

#### Removed

- Printing of the feature limit after set up of the classification. 
  
#### Known Issues

- No support currently for MuData objects. All I need to do is add a wrapper to convert a mudata object to an acceptable anndata object (just need to create
  a concat version, where a col in .var refers to modality, and another col corresponds to any "to_use"-type columns, joined together)
  
- Currently, there is only full support for two modalities: the .X, and a data_key corresponding to a .obsm location. There is partial support for AnnData 
  objects with multiple modalities in the .var, but this is not yet supported by the `mmc.utils.get_data` or `mmc.utils.marker` functions, and result in errors 
  when used for classification. In the future, there will be a data_key location referring to either the .obsm or to a column in the .var, allowing for a .var 
  column that specifies many (not just 2) modalities. 

- There is currently no validation of data_keys or feature names, but they cannot include `_obs`, `_gex`, or `_mod_`. 

### Earlier versions did not include a detailed change log. 
