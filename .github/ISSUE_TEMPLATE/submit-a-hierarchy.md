---
name: Submit a hierarchy
about: Submit a new hierarchy to be hosted in MMoCHi's documentation
title: ''
labels: submission
assignees: ''

---

**Name**
The name of the hierarchy. Please include a version number.
e.g. "Human T cell subsets v1"

**Short description**
e.g. "Classifier for 8 subsets of αβ T cells and monocytes."

**Paired modalities**
Please provide the modalities used, along with details on the suffix labels used in the hierarchy. 
e.g. "CITE-seq: gene expression (suffixed with `_gex`) and protein (no modifier)"

**Feature set details**
Details about what feature sets were used for high-confidence thresholding and classification.
e.g. "Aligned to GRCh38 with Gencode v24 annotation and includes antibodies from a custom universal TotalSeq-A panel of ~270 antibodies (BioLegend: 99786). Classified using all protein-coding genes and all proteins, excluding isotype controls." 

**Species**

**Sample or dataset details**
Please provide information about the sample (tissue type, preprocessing) and/or reference(s) for the dataset(s) this hierarchy was applied to.
e.g. "Monocytes and αβ T cells sorted by FACS from human PBMCs"

**Author(s)**

**Publication**
Is this hierarchy used in a publication? If so, please provide the DOI.

**Additional details**
Any additional notes, including references for any unique subsets/markers and descriptions for any validation steps performed. Please also include any tips for how to draw thresholds for any markers that might be tricky! If updating a hierarchy, please also provide a reason for updating. 

**Example code to build the hierarchy**
In the code block below, please provide the code used to build the hierarchy (i.e. the `h.add_classification` and `h.add_subset` commands). Please make sure that all markers referenced follow the conventions you listed under paired modalities.
```
h = mmc.Hierarchy() 
```
