"""Configuration file for the Sphinx documentation builder."""
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
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = "Maud"
copyright = (
    "2019, Novo Nordisk Foundation Center for Biosustainability, "
    "Technical University of Denmark"
)
author = (
    "Novo Nordisk Foundation Center for Biosustainability, "
    "Technical University of Denmark"
)

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinx_click.ext",
    "myst_parser",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# So that readthedocs will look for index.rst rather than contents.rst
master_doc = "index"

# Including numbered cross-referencing of figures
numfig = True

math_number_all = True
math_eqref_format = "eq{number}"
math_numfig = True

myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "smartquotes",
]

bibtex_bibfiles = ["bibliography.bib"]
bibtex_reference_style = "author_year"
bibtex_default_style = "plain"

html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
