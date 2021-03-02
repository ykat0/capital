# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import sys
from datetime import datetime
import sphinx_rtd_theme
from pathlib import Path

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent)]

import capital

# -- Project information -----------------------------------------------------

project = 'CAPITAL'
author = capital.__author__
#copyright = f'{datetime.now():%Y}, {author}'
copyright = f'{datetime.now():%Y}'
release = '1.0.0'
version = capital.__version__

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'nbsphinx'
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True
autodoc_typehints = "none"


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

source_suffix = ['.rst', '.md']

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
