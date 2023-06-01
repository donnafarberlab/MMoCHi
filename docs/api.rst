API
***

Import MMoCHi as:

```import mmochi as mmc```


Hierarchy Generation
--------------------

Designing your hierarchy is the first step to classification. The Hierarchy class contains many methods used for building the hierarchy and defining high confidence thresholds. For usage examples, see the :ref:`tutorials`.

.. module:: mmochi.hierarchy

.. autosummary::
   :toctree: functions/

   Hierarchy
   Subset
   Classification
   hc_defs
   ...

Classification 
--------------

.. module:: mmochi.classifier

.. autosummary::
   :toctree: functions/

   classifier_setup
   hc_threshold
   classify
   terminal_names
   ...

Thresholding
------------

.. module:: mmochi.thresholding

.. autosummary::
   :toctree: functions/

   threshold
   run_threshold
   ...

Landmark Resgistration 
----------------------

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

Plotting 
--------

Once you have run your classification, you may be interested in plotting some metrics of its performance, or evaluating feature importances. For that, you can use these plotting functions here. For usage examples, see the :ref:`tutorials`.

.. module:: mmochi.plotting

.. autosummary::
   :toctree: functions/

   plot_confusion
   plot_confidence
   feature_importances
   plot_important_features
   plot_tree
   ...

There are also a few helper functions we have created for interrogating high confidence thresholds and classifier performance using UMAPs: 

.. module:: mmochi.utils
   :noindex:

.. autosummary::
   :toctree: functions/

   umap_thresh
   umap_interrogate_level
   ...


Helper functions
----------------

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

.. module:: mmochi.classifier
   :noindex:

.. autosummary::
   :toctree: functions/

   identify_group_markers
   ...

Logging
----------------

.. module:: mmochi.logger

.. autosummary::
   :toctree: functions/

   log_to_file
   ...

.. tocTree::
   :hidden:
   :maxdepth: 2
   :titlesonly:



