# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Ensure src/ is on sys.path for autodoc
sys.path.insert(0, os.path.abspath("../src"))

project = "MANTIS"
copyright = "2026, Uli Wortmann"
author = "Uli Wortmann"

# Retrieve version dynamically from package if installed or fallback
try:
    from mantis._version import version as release
except Exception:
    release = "0.1.0"
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "myst_parser",
]

# MyST Markdown parser settings
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "deflist",
    "fieldlist",
    "html_image",
]
myst_heading_anchors = 3
suppress_warnings = ["myst.xref_missing"]

# Napoleon settings for NumPy-style docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_param = True
napoleon_use_rtype = True

# Autosummary settings
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"

# Pygments custom lexer mapping
from pygments.lexers import PythonLexer
from sphinx.highlighting import lexers
lexers["jupyter-python"] = PythonLexer()

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "_minted", "tutorial.md", "Tutorial_new.md", "Tutorial_new.tex", "Tutorial_new.pdf", "export_tutorials.sh"]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = "_static/mantis.svg"
html_css_files = ["custom.css"]
html_theme_options = {
    "logo_only": True,
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
}
