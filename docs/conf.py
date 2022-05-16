from unittest.mock import MagicMock

import sys
from pathlib import Path
from datetime import datetime


HERE = Path(__file__).parent


# Can’t use autodoc_mock_imports as we import anndata2ri
sys.modules["rpy2.rinterface_lib"] = MagicMock()
submods = ["embedded", "conversion", "memorymanagement", "sexp", "bufferprotocol", "callbacks", "_rinterface_capi"]
sys.modules.update({f"rpy2.rinterface_lib.{sub}": MagicMock() for sub in submods})
sexp = sys.modules["rpy2.rinterface_lib.sexp"]
sexp.Sexp = type("Sexp", (MagicMock,), dict(__module__="rpy2.rinterface"))
sexp.SexpVector = sexp.SexpEnvironment = sexp.StrSexpVector = MagicMock
sexp.SexpVector.from_iterable = MagicMock()

import rpy2.rinterface

rpy2.rinterface._MissingArgType = object

# now we can import it!
sys.path[:0] = [str(HERE.parent), str(HERE / "ext")]
import anndata2ri.scipy2ri  # noqa


# -- General configuration ------------------------------------------------


needs_sphinx = "1.7"  # autosummary bugfix

# General information
project = "anndata2ri"
author = anndata2ri.__author__
copyright = f"{datetime.now():%Y}, {author}."
version = anndata2ri.__version__.replace(".dirty", "")
release = version

# default settings
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"
# default_role = '?'
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "sphinx"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx_autodoc_typehints",
    "scanpydoc",
    *[p.stem for p in (HERE / "ext").glob("*.py")],
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
todo_include_todos = False

intersphinx_mapping = dict(
    anndata=("https://anndata.readthedocs.io/en/latest/", None),
    numpy=("https://docs.scipy.org/doc/numpy/", None),
    pandas=("http://pandas.pydata.org/pandas-docs/stable/", None),
    python=("https://docs.python.org/3", None),
    rpy2=("https://rpy2.github.io/doc/latest/html/", None),
    scipy=("https://docs.scipy.org/doc/scipy/reference/", None),
)


# -- Options for HTML output ----------------------------------------------


html_theme = "scanpydoc"
html_theme_options = dict(collapse_navigation=True)
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_context = dict(
    display_github=True,
    github_user="theislab",
    github_repo="anndata2ri",
    github_version="master",
    conf_py_path="/docs/",
)


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (master_doc, project, doc_title, author, project, "One line description of project.", "Miscellaneous")
]
