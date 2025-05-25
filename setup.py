from setuptools import setup, find_packages

VERSION = '0.3.5dev' 
DESCRIPTION = 'MultiMOdal Classifier HIerarchy (MMoCHi)'
LONG_DESCRIPTION = 'A reference-free hierarchical classification system designed for high-confidence thresholding of cell populations on CITE-Seq and training of random forests.'

setup(  name="mmochi", 
        version=VERSION,
        author="Daniel Caron",
        author_email="<dpc2136@cumc.columbia.edu>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        url = 'https://github.com/donnafarberlab/MMoCHi.git',
        keywords=['scRNA','CITE-Seq','Classifier','Multimodal'],
	    install_requires=['pandas>=1.0.5','numpy>=1.21.0','scikit-learn>=1.1.2','imbalanced-learn>=0.11.0',
                          'scipy>=1.8.0','scanpy>=1.8.0','anndata>=0.8.0','treelib>=1.6.1',
                          'matplotlib>=3.6.1','ipywidgets>=7.7.0','pydot>=1.4.2','ipython>=8.4.0',
                          'tqdm>=4.64.0','leidenalg>=0.8.9'],
        extras_require= {'landmarking':['scikit-fda>=0.5','seaborn>=0.11.2'],
                         'tutorials':['harmonypy'],
                         'docs':['sphinx==6.2.1','sphinx-rtd-theme==1.2.1','m2r2==0.3.2',
                                 'nbsphinx==0.9.2','jinja2==3.0.3','docutils','sphinx-autodoc-typehints==1.21.8',
                                 'pandoc','sphinx_new_tab_link','jupyter_sphinx','sphinx-copybutton','lxml_html_clean'],
                         'pytest':['pytest','pytest-cov','pytest-fixtures','tox','nbmake','tox-conda']},
        classifiers= [])