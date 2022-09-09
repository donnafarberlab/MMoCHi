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

MMoCHi requires Python == 3.8

The easiest way to setup an environment is with Anaconda Python distribution in Miniconda or anaconda:
```
conda create -n mmochi python=3.8
```
To activate this environment in jupyter, you can either activate the conda environment before starting jupyter:
```
conda activate mmochi
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
Note that with this option, if you are in jupyter labs <3.4.2, running terminal commands (such as `pip install`) via the notebook will still open a new terminal in the base environment, so editing the enviornment will require opening an new terminal window and running `conda activate mmochi` to enter the environment.

To use the built-in landmark registration functions, you will also need to install scikit-fda:
```
conda activate mmochi
pip install scikit-fda==0.5
```


While not required, the following dependencies can be installed to enable the use of plotting functions for displaying complex hierarchies or trees from a random forest:

```
conda activate mmochi
pip install pydot
sudo apt-get install graphviz
```



To follow along in our tutorial, we also recommend installing the beta version of pySankey:

```
conda activate mmochi
pip install pySankeyBeta
```

Lastly, future versions of MMoCHi will support MuData objects, which can be installed with the following:

```
conda activate mmochi
pip install mudata
```

### Installing from source
Once you have set up the environment, clone this repository and install:
```
git clone https://github.com/donnafarberlab/mmochi.git
cd mmochi
pip install .
```

### Testing your installation
This will verify a successful install of mmochi. Testing will catch major issues in dependencies. From the mmochi base
directory do:
```
conda install pytest
pytest
```
<div class="alert alert-block alert-warning">
<b>Testing Incomplete</b> Testing is currently only set up for a few basic functions, more robust testing is on the way! For now, running through the Classifier Demo tutorial will help ensure your environment is correct.
</div>

## Quick Start: Tutorials

- [Integrated Classification](/docs/Classifier_Demo.ipynb) - Learn how to set up and run MMoCHi on 10X Genomics CITE-Seq!

<div class="alert alert-block alert-warning">
<b>Under Construction</b> More tutorials are on the way, both more basic and more advanced. Eventually hosting on collab may also become possible.
</div>

## Docs

<div class="alert alert-block alert-danger">
<b>In progress</b> Robust docs are on the way, but for now, most functions' docstrings are well documented. Type hints are also in progress. 
</div>

## Known issues

- `_in_danger_noise` function is currently computationally unoptimized
- Not all dependencies are necessary, and some need to be made conditional upon function use.
- No support currently for MuData objects, or other ways to store multimodal data.
- Log files are currently saved and overwrite immediately upon package loading, with no option to change logging or verbosity.
