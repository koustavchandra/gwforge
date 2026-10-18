# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "GWForge"
copyright = "2023, Koustav Chandra"
author = "Koustav Chandra"
# Single source of truth is setup.cfg; the Pages build installs the package first.
from importlib.metadata import version as _version

release = _version("gwforge")

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # myst_nb renders the tutorial notebooks and bundles myst_parser for the pages.
    "myst_nb",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    # Draws the dependency flowchart on the code-overview page.
    "sphinxcontrib.mermaid",
]
extensions.append("autoapi.extension")
autoapi_dirs = ["../../GWForge"]
# Render numpydoc Attributes as :ivar: fields; autoapi already documents the
# attributes themselves, and the two together are "duplicate object" warnings.
napoleon_use_ivar = True
# Document each object where it is defined, as sphinx-apidoc would; the
# package pages keep their summary tables and link down.
# The API entry is placed by hand under its own caption in index.md.
autoapi_add_toctree_entry = False
autoapi_options = ["members", "undoc-members", "show-inheritance", "show-module-summary"]
# The pages are MyST Markdown: `$...$` and `$$...$$` are parsed as maths here,
# before the Markdown parser can touch the LaTeX inside them.
myst_enable_extensions = ["dollarmath", "amsmath"]
myst_heading_anchors = 3
# The tutorials are executed locally and committed with their outputs: they
# need lalsuite and minutes of runtime, which the Pages build does not have.
nb_execution_mode = "off"
# Logging and leap-second warnings are stderr noise on a rendered page.
nb_output_stderr = "remove"
autosectionlabel_prefix_document = True
exclude_patterns = ["**.ipynb_checkpoints"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
