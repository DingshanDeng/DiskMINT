# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# If your package uses a src/ layout: repo/src/diskmint
# conf.py is at repo/docs/source/conf.py, so go up twice to reach repo/
# (optional; not required for autodoc2, but helpful if you import anywhere)
sys.path.insert(0, os.path.abspath(os.path.join(__file__, "..", "..", "..", "diskmint", "src")))
PKG_PATH = os.path.abspath(os.path.join(__file__, "..", "..", "..", "diskmint", "src", "diskmint"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DiskMINT'
copyright = '2023-2025, Dingshan Deng'
author = 'Dingshan Deng'
# release = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [
#   'myst_parser', # MyST Markdown
#   'sphinx_design',
#   'autodoc2',
# #   'sphinx.ext.autodoc',
#   'sphinx.ext.coverage',
#   'sphinx.ext.napoleon',
#   'sphinx.ext.imgmath',
#   'nbsphinx',
#   'IPython.sphinxext.ipython_console_highlighting',
#   "sphinx.ext.mathjax",   # MathJax support
# ]

extensions = [
    'myst_parser',
    'sphinx_design',
    'autodoc2',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',        # prefer MathJax over imgmath on RTD
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
]

master_doc = 'index'

# Get that autodoc working.
# autodoc_mock_imports = ['astropy', 'scipy', 'numpy', 'pandas', 'matplotlib']
# autodoc_member_order = 'bysource'

# autodoc2_packages = ['../../src/diskmint']

# autodoc2: tell it where your package lives and to parse Markdown docstrings
autodoc2_packages = [
    {"path": PKG_PATH, "module": "diskmint", "auto_mode": True},
]

# Put the generated API pages under docs/apidocs/
autodoc2_output_dir = "apidocs"

# Parse & render docstrings as Markdown via MyST
autodoc2_docstring_parser = "myst"
autodoc2_render_plugin = "myst"

# Useful filters (optional)
autodoc2_hidden_objects = ["dunder", "private"]   # hide __dunder__ and _private
# autodoc2_module_all = True                      # respect __all__ if you use it

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

myst_enable_extensions = [
    'amsmath',  # LaTeX math environments
    'attrs_inline',
    'colon_fence',
    'deflist',
    'dollarmath', # inline $...$ and display $$...$$
    'fieldlist',
    'html_admonition',
    'html_image',
    'linkify',
    'replacements',
    'smartquotes',
    'strikethrough',
    'substitution',
    'tasklist',
    # "attrs",  
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints', '.ideas']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# Readthedocs.
# on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
# if not on_rtd:
#     import sphinx_rtd_theme
#     html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named 'default.css' will overwrite the builtin 'default.css'.
html_static_path = ['_static']
html_css_files = ["css/custom.css"]   # make this file in docs/_static/css/custom.css

# html_logo = '_static/logo_mini.png'
html_logo = '_static/assets/images/card-software-transparent.png'
html_theme_options = {
    "logo_only": True,       # either is fine; does not block the ::before/::after text
    "display_version": False, # we’re adding our own caption via CSS anyway
    "navigation_depth": 3, # make the sidebar expand more levels
}
