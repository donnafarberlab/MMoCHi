# MMoCHi: <ins>M</ins>ulti<ins>Mo</ins>dal <ins>C</ins>lassifier <ins>Hi</ins>erarchy

<img align="right" src="./docs/_static/mmochi_logo.svg" width="480">

## About
MMoCHi is a hierarchical classification system designed around CITE-Seq data, but applicable to diverse single cell sequencing modalities. 
It includes tools for:

1. <b>Integrating data</b> of Antibody Derived Tag (ADT) expression
2. <b>High confidence thresholding</b> of cell populations on ADT expression, gene expression, and cell metadata 
3. <b>Building a hierarchy</b> of cell subset classifiers 
4. <b>Training and evaluating</b> random forest classifiers at each hierarchy level

## Docs
A ReadTheDocs is on the way! For now, refer to docstrings of each function.

## Tutorials

### Quick Start

- **[Integrated Classification](/docs/Integrated_Classification.ipynb)** - Set up and run MMoCHi on 10X Genomics CITE-Seq data!
- **[Input/Output Specifications](/docs/Input_Output_Specs.md)** - A handy reference for expected input formats and the meaning of various outputs.

### In-Depth Resources 

- **[Landmark Registration](./docs/Landmark_Registration.ipynb)** - Improving integration by manually adjusting peak identification.
- **[Hierarchy Design](./docs/Hierarchy_Design.ipynb)** - Building your own MMoCHi hierarchy object and selecting markers for high confidence definitions.
- **[High Confidence Thresholding](./docs/High_Confidence_Thresholding.ipynb)** - Careful thresholding to improve high confidence selection and training!
- **[Exploring Feature Importances](./docs/Exploring_Feature_Importances.ipynb)** - Interpreting random forests: What features are most useful for classification?
- **[Pretrained Classification](./docs/Pretrained_Classification.ipynb)** - Applying a pretrained classifier to other datasets.

# Installation Instructions

### Setting up a Conda Environment

MMoCHi is built for Python == 3.8, but may also work on later python versions.

The easiest way to setup an environment is with a Python distribution in miniconda or anaconda:
```
conda create -n mmochi python=3.8
```
### Options to set up Jupyter

MMoCHi was designed to be run in [Jupyter IPython Notebooks](https://jupyter.org/). To use this environment in Jupyter Labs you can either create a fresh install of Jupyter within this new conda environment:
```
conda activate mmochi
conda install -c conda-forge jupyterlab
jupyter lab
```
or you can install this new conda environment as a new kernel to your current Jupyter installation.
```
conda activate mmochi
conda install ipykernel
python -m ipykernel install --user --name=mmochi
conda deactivate
juptyer lab
```
Note that with this option, the mmochi environment is accessed via the "Kernel>Change kernel..." menu. Note that with this option, depending on your Jupyter version, terminal commands (such as `pip install`) run within a code block may be run in the base environment by default, so editing the conda environment will require running `conda activate mmochi` in a new terminal window to enter the environment.

### Cloning the Repository
Once you have set up the environment, clone this repository and `cd` into its root directory.:
```
conda activate mmochi
git clone https://github.com/donnafarberlab/mmochi.git
cd mmochi
```
Alternatively, you can download the repository as a .zip file, extract it, and navigate to its root directory (likely `MMoCHi-main`).

### Installing Dependencies

Most of the dependencies for MMoCHi can be installed via conda using:
```
conda activate mmochi
conda install -c conda-forge --file requirements.txt
```

While not required, graphviz can be installed to enable the use of plotting functions for displaying complex hierarchies or trees from a random forest, and harmonypy can be installed to follow along in our tutorial:
```
conda activate mmochi
sudo apt-get install graphviz
pip install harmonypy
```

If you are on an Intel(R) processor, you can also install the following package to accelerate computation:
```
conda activate mmochi
pip install scikit-learn-intelex
```

### Installing MMoCHi
Lastly, install the package from source:
```
conda activate mmochi
pip install .
```

Optional dependencies for landmark registration, compiling documentation, or running pytest can be installed using one of the following:
```
pip install .[landmarking]
pip install .[docs]
pip install .[pytest]
```

### Compiling the docs [OPTIONAL]

After running `pip install .[docs]`, enter the docs directory and compile the docs:

```
cd docs
make html
```

There will be many warnings, but it should end with `"build succeeded"`. You can then view the docs by opening 
`docs/_build/html/index.html`, which should pop up as an interactive site in your local browser.

If you are having trouble with your pandoc installation, try uninstalling it with pip, and reinstalling it with conda.

```
pip uninstall pandoc
conda install -c conda-forge pandoc
```


### Testing your installation
This will help verify a successful install of MMoCHi. While the pytest covers ~90% of the code, it still does not capture many use-cases. Thus, while it will catch major issues in dependencies, it may not detect issues with functions requiring user interaction. Please also run through the tutorials to verify your environment is correct.
From the mmochi folder run:
```
conda activate mmochi
conda install pytest
pytest
```
