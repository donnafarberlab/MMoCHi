from setuptools import setup, find_packages

VERSION = '0.1.3' 
DESCRIPTION = 'MultiMOdal Classifier HIerarchy (MMoCHi)'
LONG_DESCRIPTION = 'A hierarchical classification system designed for ground-truth thresholding of cell populations on CITE-Seq for training of random forests.'

setup(  name="mmochi", 
        version=VERSION,
        author="Daniel Caron",
        author_email="<dpc2136@cumc.columbia.edu>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        url = 'https://github.com/donnafarberlab/MMoCHi.git',
        keywords=['scRNA','CITE-Seq','Classifier','Multimodal'],
	    install_requires=['pandas','numpy','scikit-learn','imbalanced-learn',
                          'scipy','scanpy','anndata','treelib','matplotlib'],
        classifiers= [])
