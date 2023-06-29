# MMoCHi: <ins>M</ins>ulti<ins>Mo</ins>dal <ins>C</ins>lassifier <ins>Hi</ins>erarchy

<img align="right" src="./docs/_static/mmochi_logo.svg" width="480">

## About
MMoCHi is a hierarchical classification system designed for CITE-Seq data, but applicable to diverse single-cell sequencing modalities. 
It includes tools for:

1. <b>Integrating data</b> of Antibody Derived Tag (ADT) expression
2. <b>High-confidence thresholding</b> of cell populations using ADT expression, gene expression, and cell metadata 
3. <b>Building a hierarchy</b> of cell subset classifiers 
4. <b>Training and evaluating</b> random forest classifiers at each hierarchy level

## Docs
Documentation is on the way!

## Tutorials

### Quick Start

- **[Integrated Classification](/docs/Integrated_Classification.ipynb)** - Set up and run MMoCHi on 10X Genomics CITE-Seq data!
- **[Input/Output Specifications](/docs/Input_Output_Specs.md)** - A handy reference for expected input formats and the meaning of various outputs.

### In-Depth Resources 

- **[Landmark Registration](./docs/Landmark_Registration.ipynb)** - Improving integration by manually adjusting peak identification.
- **[Hierarchy Design](./docs/Hierarchy_Design.ipynb)** - Building your own MMoCHi hierarchy object and selecting markers for high-confidence definitions.
- **[High-Confidence Thresholding](./docs/High_Confidence_Thresholding.ipynb)** - Careful thresholding to improve high-confidence selection and training!
- **[Exploring Feature Importances](./docs/Exploring_Feature_Importances.ipynb)** - Interpreting random forests: What features are most useful for classification?
- **[Pretrained Classification](./docs/Pretrained_Classification.ipynb)** - Applying a pretrained classifier to other datasets.

# Installation Instructions

### Setting up a virtual environment

MMoCHi was developed using Python 3.8, but may also work on later Python versions. We recommend isolating installation using a virtual environment. You can create a virtual environment using miniconda or anaconda:
```
conda create -n mmochi python=3.8
```
Once you create this enviornment, activate it: 
```
conda activate mmochi
```
Ensure this environment is activated for all of the following installation steps. For this reason, we recommend running installation in terminal, not through IPython notebooks. 

### Setting up Jupyter

MMoCHi was designed to be run in [Jupyter IPython Notebooks](https://jupyter.org/). If you already have JupyterLab installed, you can use ipykernel to install your new environment into Jupyter as a separate kernel. 
```
conda install ipykernel
python -m ipykernel install --user --name=mmochi
```
This environment can now be accessed within Jupyter via the "Kernel>Change kernel..." menu. You may need to restart Jupyter for changes to apply.

**ALTERNATIVELY**, if you do not have Jupyter preinstalled, you could create a fresh install of JupyterLab within this new conda environment:
```
conda install -c conda-forge jupyterlab
jupyter lab
```

### Cloning the repository

Once you have set up the environment, clone this repository and `cd` into its root directory: 
```
git clone https://github.com/donnafarberlab/mmochi.git
cd mmochi
```
(Note: The root directory is the one the with README.md in it)

**ALTERNATIVELY**, you can download the repository as a .zip file, extract it, and navigate to its root directory (which may be named `MMoCHi-main`).

### Installing dependencies
Most of the dependencies for MMoCHi can be installed via conda:
```
conda install -c conda-forge --file requirements.txt
```

While not required, additional dependencies can be installed to enable hierarchy plotting functions (graphviz) and to follow along in our tutorials (harmonypy):
```
sudo apt-get install graphviz
python -m pip install harmonypy
```

If you are on an Intel(R) processor, you can also install the following package to accelerate computation:
```
python -m pip install scikit-learn-intelex
```

### Installing MMoCHi
Lastly, while in MMoCHi's root directory, install the package from source. This install also includes dependencies for running landmark registration and for testing your installation with pytest.
```
python -m pip install ".[landmarking,pytest]"
```

### Testing your installation
Testing can help verify a successful install of MMoCHi. While the pytest covers ~90% of the code, it still does not capture many use-cases. Thus, while it catches major dependency issues, it may not detect issues with functions requiring user interaction. Please also run through the tutorials to verify your environment is correct. From MMoCHi's root directory run:
```
conda activate mmochi
python -m pytest
```

### Troubleshooting installation
Although conda usually handles conflicts between dependency versions, issues sometimes still arise. Here, we have also provided a list of set dependencies for Python 3.8. Recreate your environment, set up Jupyter, navigate to MMoCHi's root directory, then install as follows:
```
conda activate mmochi
python -m pip install -r Python3_8_requirements.txt
python -m pip install .
```

### For developers
Contributions to MMoCHi are welcome. Please get in touch if you would like to discuss. To install for development, navigate to MMoCHi's root directory and create an editable install as follows:
```
python -m pip install -e ".[docs,landmarking,pytest,tutorials]"
```
In addition to running the testing suite, you can also test tutorial notebooks automatically, check code coverage, and test many python versions and dependency sets:
```
python -m pytest --nbmake docs/Integrated_Classification.ipynb docs/Landmark_Registration.ipynb docs/Hierarchy_Design.ipynb docs/High_Confidence_Thresholding.ipynb docs/Exploring_Feature_Importances.ipynb docs/Pretrained_Classification.ipynb

python -m pytest --cov=mmochi --cov-report=html 
python -m pytest --nbmake docs/*.ipynb --cov=mmochi --cov-report=html --cov-append

python -m tox -r
```

You can also compile the docs using the following command:
```
cd docs
make html
```
There will be many warnings, but it should end with `"build succeeded"`. You can then view the docs by opening `docs/_build/html/index.html` in your local browser.

Occationally, there are errors running pandoc when installed using pip. Try resinstalling with conda:
```
pip uninstall pandoc
conda install -c conda-forge pandoc
cd docs
make html
```



