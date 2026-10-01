# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "GWForge"
copyright = "2023, Koustav Chandra"
author = "Koustav Chandra"
release = "0.0.dev1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
]
extensions.append("autoapi.extension")
autoapi_dirs = ["../../GWForge"]
# The pages are MyST Markdown: `$...$` and `$$...$$` are parsed as maths here,
# before the Markdown parser can touch the LaTeX inside them.
myst_enable_extensions = ["dollarmath", "amsmath"]
myst_heading_anchors = 3
autosectionlabel_prefix_document = True
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
