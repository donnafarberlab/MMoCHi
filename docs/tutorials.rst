Tutorials
*********


:ref:`**Integrated Classification with MMoCHi**`

After :ref:`installation`, get started learning how to use MMoCHi by following our integrated classification tutorial for an overview integration, hierarchy creation, classifier training, and classifier application on a CITE-Seq data.


:ref:`**Landmark Registration — using the built in MMoCHi function**`

Landmark registration is a useful inhouse tool for integration of ADT expression across batches. This tutorial details how to location positive and negative populations for integration and how to tune landmark registration to obtain well-integrated ADT expression. 


:ref:`**Hierarchy Design for MMoCHi**`

Careful hierarchy design is a critical part of effective classification. This tutorial explains different ways to define a MMoCHi hierarchy and best practices to achieve effective annotation.


:ref:`**High-Confidence Thresholding**`

High-confidence thresholding is the process by which MMoCHi selects training data for its classifier. MMoCHi uses expression thresholding, an approach similar to flow cytometry, on gene and protein expression to select events for training. This tutorial walks through validation and tuning of those thresholds.


:ref:`**Exploring Feature Importances with MMoCHi**`

Once a classifier has been trained, MMoCHi offers insight into the decisions made at each level of classification by exploring feature importances for performance insights and identification of key gene and protein markers that differentiate cell types. This tutorial walks through using these functions and some of the conclusions that can be drawn from certain important features. 


:ref:`**Pretrained Classification using MMoCHi**`

MMoCHi can also be used for label transfer using a pretrained classification. The pretrained classification tutorial demonstrates how to take a classifier and apply it directly to a new dataset to achieve this goal.

:ref:`**Automatic Hyperparameter Optimization with MMoCHi**`

Random forests in MMoCHi can sometimes require hyperparameter tuning. This can be done manually, or using MMoCHi's built-in functionality.

.. tocTree::
   :maxdepth: 0
   :titlesonly:
   :hidden:

   Integrated_Classification
   Landmark_Registration
   Hierarchy_Design
   High_Confidence_Thresholding
   Exploring_Feature_Importances
   Pretrained_Classification
   Hyperparameter_Optimization