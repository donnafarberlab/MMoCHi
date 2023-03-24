# MMoCHi: <ins>M</ins>ulti<ins>Mo</ins>dal <ins>C</ins>lassifier <ins>Hi</ins>erarchy

<img align="right" src="./docs/figures/mmochi_logo.svg" width="480">

## About

MMoCHi is a <i>hierarchical classification</i> system designed around CITE-Seq data, but applicable to diverse single cell sequencing modalities. 
It includes tools for:

1. <b>Integrating data</b> of Antibody Derived Tag (ADT) expression
2. <b>Ground-truth thresholding</b> of cell populations on ADT expression, gene expression, as well as cell metadata 
3. <b>Building a hierarchy</b> of cell subset classifiers 
4. <b>Training and evaluating</b> random forests at each hierarchy level

## Installation

### Environment and Dependencies

MMoCHi is built for Python == 3.8, but may also work on later python versions.

The easiest way to setup an environment is with Anaconda Python distribution in Miniconda or anaconda:
```
conda create -n mmochi python=3.8
```
To activate this environment in jupyter, you can either activate the conda environment before starting jupyter:
```
conda activate mmochi
conda install -c conda-forge jupyterlab
jupyter lab
```
or you can install your conda environment as an ipykernel, then access it via the "Kernel>Change kernel..." menu
```
conda activate mmochi
conda install ipykernel
python -m ipykernel install --user --name=mmochi
conda deactivate
juptyer lab
```
Note that with this option, if you are in jupyter labs <3.4.2, running terminal commands (such as `pip install`) via the notebook will still open a new terminal
in the base environment, so editing the enviornment will require opening a new terminal window and running `conda activate mmochi` to enter the environment.

To use the built-in landmark registration functions, you will also need to install scikit-fda:
```
conda activate mmochi
pip install scikit-fda==0.5
```

While not required, the following dependencies can be installed to enable the use of plotting functions for displaying complex hierarchies or trees from a 
random forest:

```
conda activate mmochi
pip install pydot
sudo apt-get install graphviz
```

To follow along in our tutorial and use all features, we also recommend installing harmonypy, leidenalg, the beta version of pySankey, and iPyWidgets:

```
conda activate mmochi
pip install harmonypy
pip install leidenalg
pip install pySankeyBeta
pip install ipywidgets
```

### Installing from source
Once you have set up the environment, clone this repository and install:
```
conda activate mmochi
git clone https://github.com/donnafarberlab/mmochi.git
cd mmochi
pip install .
```

### Testing your installation
This will verify a successful install of mmochi. Testing will catch major issues in dependencies. From the mmochi base directory do:
```
conda activate mmochi
conda install pytest
pytest
```
<div class="alert alert-block alert-warning">
<b>Testing Incomplete</b> While the pytest covers ~90% of the code, it still does not capture many use-cases. 
    Please also run through the demo tutorial to verify your environment is correct.
</div>

## Quick Start: Tutorials
Currently, tutorials are available as iPython notebooks. Eventually hosting on collab may also become possible.

- [Input/Output Specifications](/docs/Input_Output_Specs.md) - Detailed information about the expected inputs for running MMoCHi!

- [Integrated Classification](/docs/Classifier_Demo.ipynb) - Learn how to set up and run MMoCHi on 10X Genomics CITE-Seq!

- [Fine-tune Landmark Registration](/docs/Landmark_Registration.ipynb) - Manually adjusting the peak identification to improve integration!

- Improving Ground Truth (**In progress**) - Carefully picking markers and thresholding them can improve ground truth selection!

- Explore Feature Importances (**In progress**) - What features are most useful for classification?

- Troubleshoot Classification (**In progress**) - What can you do to improve classification performance?

## Docs

<div class="alert alert-block alert-danger">
<b>In progress</b> A ReadTheDocs or alternative is on the way, but for now, most functions have clear docstrings and type hints.
</div>