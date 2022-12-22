# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/) and after v1.0.0 release, 
this project will adhere to [Semantic Versioning](http://semver.org/).
 
# Current version
## [0.1.3] - 22Dec22

To update, replace your repository with a fresh clone, and install using pip, as before.

```
conda activate mmochi
git clone https://github.com/donnafarberlab/mmochi.git
cd mmochi
pip install .
```

### Added
 
- Updated tutorials for [Integrated Classification](/docs/Classifier_Demo.ipynb) and [Fine-tuning Landmark Registration](/docs/Landmark_Registration.ipynb)
  
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

### Changed
  
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
 
### Fixed

- Improved pytest coverage of all modules, to test broader use-cases.

- Specified peak overrides in `mmc.landmarking._landmark_registration` should now more accurately reflect the closest possible mapping of those values to the 
  proper position in the FDA function.

- Using `return_graph=True` and `plot=False` together on `mmc.hierarchy.Hierarchy.display` will now return None (expected behavior) instead of throwing an
  error. 

- When initializing a new `mmc.hierarchy.Classification` with a predefined classifier, there is now error checking to ensure feature_names are also provided

- Warnings about any cython import failures during the import of scikit-fda are now silenced.

- Fixed performance of `mmc.utils.generate_exlcusive_features` if `adatas` is given as a list of str. Previously this would mistakenly return an empty list.

- Fixed reading in 10x formatted .h5 files without url backend using `mmc.utils.preprocess_adatas`

- Fixed `get_data` function when passed a keyword including `"_obs"` in the variable name

- Added fixes for cases where thresholds in `umap_thresh` were out of the bounds of the data.

- Fixed error where if features were limited on a per-classification-level basis, the wrong set of features were passed to generate the training matricies.

- Fixed error if no clusters needed to be balanced in `mmc.classifier.balance_borderline`

### Removed

- Printing of the feature limit after set up of the classification. 
  
### Known Issues

- No support currently for MuData objects. All I need to do is add a wrapper to convert a mudata object to an acceptable anndata object (just need to create
  a concat version, where a col in .var refers to modality, and another col corresponds to any "to_use"-type columns, joined together)
  
- Currently, there is only full support for two modalities: the .X, and a data_key corresponding to a .obsm location. There is partial support for AnnData 
  objects with multiple modalities in the .var, but this is not yet supported by the `mmc.utils.get_data` or `mmc.utils.marker` functions, and result in errors 
  when used for classification. In the future, there will be a data_key location referring to either the .obsm or to a column in the .var, allowing for a .var 
  column that specifies many (not just 2) modalities. 

- There is currently no validation of data_keys or feature names, but they cannot include `_obs`, `_gex`, or `_mod_`. 

# Past versions

Past versions did not include a detailed change log. 