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
import sys
import inspect
from pathlib import Path
from ipywidgets.embed import DEFAULT_EMBED_REQUIREJS_URL

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent)]
import mmochi

html_js_files = [
    DEFAULT_EMBED_REQUIREJS_URL,
]

project = 'MMoCHi'
copyright = '2023, Daniel Caron'
author = 'Daniel Caron'
repository_url = "https://github.com/donnafarberlab/MMoCHi"

version = mmochi.__version__
release = version

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = ['sphinx.ext.napoleon','m2r2','sphinx.ext.autodoc','sphinx_autodoc_typehints',
              'sphinx_rtd_theme','nbsphinx','sphinx.ext.autosectionlabel','sphinx.ext.autosummary',
              'sphinx_new_tab_link', 'sphinx.ext.intersphinx','sphinx_copybutton']

source_suffix = ['.rst','.md']
templates_path = ['_templates']

nbsphinx_execute = 'auto'
autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_typehints_format = 'fully-qualified'
autodoc_preserve_defaults = True

nbsphinx_widgets_path = '' # Prevents duplicate widgets

html_context = dict(
    display_github=True,      # Integrate GitHub
    github_user='donnafarberlab',   # Username
    github_repo='MMoCHi',     # Repo name
    github_version='main',  # Version
    conf_py_path='/docs/',    # Path in the checkout to the docs root
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
html_favicon='_static/mmochi_logo.svg'

html_css_files = [
    'css/custom.css',
]


html_theme_options={'logo_only':True,'display_version':True}#,'style_nav_header_background':'#ff7b00'}



# The following code as well as base.rst for linking to sourcecode from functions via edit on Github is inspired by https://github.com/readthedocs/readthedocs-sphinx-ext/issues/28

def get_obj_module(qualname):
    """Get a module/class/attribute and its original module by qualname"""
    modname = qualname
    classname = None
    attrname = None
    while modname not in sys.modules:
        attrname = classname
        modname, classname = modname.rsplit('.', 1)

    # retrieve object and find original module name
    if classname:
        cls = getattr(sys.modules[modname], classname)
        modname = cls.__module__
        obj = getattr(cls, attrname) if attrname else cls
    else:
        obj = None

    return obj, sys.modules[modname]


def get_linenos(obj):
    """Get an object’s line numbers"""
    try:
        lines, start = inspect.getsourcelines(obj)
    except TypeError:  # obj is an attribute or None
        return None, None
    else:
        return start, start + len(lines) - 1

#pip install must be editable for this to work
project_dir = Path(__file__).parent.parent  # project/docs/conf.py/../.. → project/
github_url = 'https://github.com/{github_user}/{github_repo}/tree/{github_version}'.format_map(html_context)
def modurl(qualname):
    """Get the full GitHub URL for some object’s qualname"""
    obj, module = get_obj_module(qualname)
    path = Path(module.__file__).relative_to(project_dir)
    start, end = get_linenos(obj)
    fragment = '#L{}-L{}'.format(start, end) if start and end else ''
    return '{}/{}{}'.format(github_url, path, fragment)


# Modify default filters to edit html_context of autosummary templates
from jinja2.defaults import DEFAULT_FILTERS

DEFAULT_FILTERS['modurl'] = modurl

# link docstrings to obj apis, inspired by scanpys approach
intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/stable/', None),
    matplotlib=('https://matplotlib.org/stable/', None),
    numpy=('https://numpy.org/doc/stable/', None),
    pandas=('https://pandas.pydata.org/pandas-docs/stable/', None),
    scipy=('https://docs.scipy.org/doc/scipy/', None),
    seaborn=('https://seaborn.pydata.org/', None),
    sklearn=('https://scikit-learn.org/stable/', None),
    scanpy=('https://scanpy.readthedocs.io/en/stable/', None),
    python=('https://docs.python.org/3', None),
)



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
