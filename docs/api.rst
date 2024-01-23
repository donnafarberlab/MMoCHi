API
***

Import MMoCHi as:

```import mmochi as mmc```

For usage examples, see the :ref:`tutorials`.


Hierarchy Generation
--------------------

Designing your hierarchy is the first step to classification. The Hierarchy class contains many methods used for building the hierarchy and defining high-confidence thresholds. 

.. module:: mmochi.hierarchy

.. autosummary::
   :toctree: functions/
   :nosignatures:
   
   Hierarchy
   Classification
   Subset
   hc_defs
   ...
   
   
Thresholding
------------

High-confidence thresholding is performed primarily through methods in the `mmc.Hierarchy` object, but these are the functions that perform thresholding under the hood. 

.. module:: mmochi.thresholding

.. autosummary::
   :toctree: functions/

   threshold
   run_threshold
   ...

Classification 
--------------

Once the Hierarchy is created and high-confidence thresholds are drawn, you are ready to classify. The `mmc.classify` function runs `mmc.classifier_setup` and `mmc.hc_threshold` internally, but you have the option to run these separately for testing.

.. module:: mmochi.classifier

.. autosummary::
   :toctree: functions/

   classifier_setup
   hc_threshold
   classify
   terminal_names
   define_external_holdout
   ...

Plotting 
--------

Once you have run your classification, you may be interested in plotting some metrics of its performance or evaluating feature importances.

.. module:: mmochi.plotting

.. autosummary::
   :toctree: functions/

   plot_confusion
   plot_confidence
   feature_importances
   plot_important_features
   plot_tree
   ...

There are also a few plotting functions we have created for interrogating high-confidence thresholds and classifier performance using UMAPs: 

.. currentmodule:: mmochi.utils

.. autosummary::
   :toctree: functions/

   umap_thresh
   umap_interrogate_level
   ...

Landmark Resgistration 
----------------------

Prior to classification, you may be interested in performing batch correction on ADT expression. This module contains the tools necessary to perform and evaluate batch correction by landmark registration.

.. module:: mmochi.landmarking

.. autosummary::
   :toctree: functions/

   landmark_register_adts
   update_landmark_register
   stacked_density_plots
   density_plot
   density_plot_total
   update_peak_overrides
   save_peak_overrides
   load_peak_overrides
   ...

Helper functions
----------------

We have also developed a suite of helper functions which may be useful for running MMoCHi or preparing your data. 

.. module:: mmochi.utils

.. autosummary::
   :toctree: functions/

   marker
   get_data
   convert_10x
   obsm_to_X
   preprocess_adatas
   intersect_features
   generate_exclusive_features
   batch_iterator
   ...

.. currentmodule:: mmochi.classifier

.. autosummary::
   :toctree: functions/

   identify_group_markers
   ...

Logging
----------------

MMoCHi has built-in logs which can be helpful for debugging or reproducibility.

.. module:: mmochi.logger

.. autosummary::
   :toctree: functions/

   log_to_file
   ...

.. tocTree::
   :hidden:
   :titlesonly: