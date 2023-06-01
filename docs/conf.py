# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# import mmochi

project = 'MMoCHi'
copyright = '2023, Daniel Caron'
author = 'Daniel Caron'
repository_url = "https://github.com/donnafarberlab/MMoCHi"

# version = mmochi.__version__
# release = version
release = '0.2.0'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.napoleon','m2r2','sphinx.ext.autodoc','sphinx_autodoc_typehints','sphinx_rtd_theme','nbsphinx',
              'sphinx.ext.autosectionlabel','sphinx.ext.autosummary']

source_suffix = ['.rst','.md']
templates_path = ['_templates']

autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_typehints_format = 'fully-qualified'
autodoc_preserve_defaults = True

html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='donnafarberlab',   # Username
    github_repo='MMoCHi',     # Repo name
    github_version='master',  # Version
    conf_py_path='/',    # Path in the checkout to the docs root
)


napoleon_google_docstring = False
napoleon_numpy_docstring = True

napoleon_use_rtype = True 
napoleon_use_param = True
html_title = "mmochi"

typehints_defaults = 'braces'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '.ipynb_checkpoints']


# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo='_static/mmochi_logo.svg'
html_theme_options={'logo_only':True,'display_version':False}#,'style_nav_header_background':'#ff7b00'}



nbsphinx_prolog = """
.. raw:: html

    <style>
        h1 {
            color: orange;
        }
        
        h2 {
            color: orange;
        }
    </style>
"""