from typing import Tuple
from unittest.mock import MagicMock

from docutils import nodes
from sphinx.application import Sphinx
from sphinx.environment import BuildEnvironment
from sphinx.roles import XRefRole


import sys
from pathlib import Path
from datetime import datetime


HERE = Path(__file__).parent


# Can’t use autodoc_mock_imports as we import anndata2ri
sys.modules["rpy2.rinterface_lib"] = MagicMock()
submods = ["embedded", "conversion", "memorymanagement", "sexp", "bufferprotocol", "callbacks", "_rinterface_capi"]
sys.modules.update({f"rpy2.rinterface_lib.{sub}": MagicMock() for sub in submods})
sexp = sys.modules["rpy2.rinterface_lib.sexp"]
sexp.Sexp = sexp.SexpVector = sexp.SexpEnvironment = sexp.StrSexpVector = MagicMock
sexp.SexpVector.from_iterable = MagicMock()

import rpy2.rinterface

rpy2.rinterface._MissingArgType = object

# now we can import it!
sys.path.insert(0, str(HERE.parent))
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


html_theme = "sphinx_rtd_theme"
html_theme_options = dict(collapse_navigation=True)
html_context = dict(
    display_github=True,
    github_user="theislab",
    github_repo="anndata2ri",
    github_version="master",
    conf_py_path="/docs/",
)


# -- Add R links ----------------------------------------------------------


class RManRefRole(XRefRole):
    nodeclass = nodes.reference

    def __init__(self, *a, cls: bool = False, **kw):
        super().__init__(*a, **kw)
        self.cls = cls

    def process_link(
        self, env: BuildEnvironment, refnode: nodes.reference, has_explicit_title: bool, title: str, target: str
    ) -> Tuple[str, str]:
        qualified = not target.startswith("~")
        if not qualified:
            target = target[1:]
        package, symbol = target.split("::")
        title = target if qualified else symbol
        topic = symbol
        if self.cls:
            topic += "-class"
        target = f"https://www.rdocumentation.org/packages/{package}/topics/{topic}"
        refnode["refuri"] = target
        return title, target

    # def result_nodes(self, document: nodes.document, env: BuildEnvironment, node: nodes.reference, is_ref: bool):
    #    target = node.get('reftarget')
    #    if target:
    #        node.attributes['refuri'] = target
    #    return [node], []


def setup(app: Sphinx):
    app.add_role("rman", RManRefRole())
    app.add_role("rcls", RManRefRole(cls=True))


# -- Quick fixes ---------------------------------------------------------------


from rpy2.rinterface import Sexp

Sexp.__module__ = "rpy2.rinterface"


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (master_doc, project, doc_title, author, project, "One line description of project.", "Miscellaneous")
]
