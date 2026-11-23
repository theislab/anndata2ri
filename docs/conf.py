"""Sphinx configuration."""

from __future__ import annotations

import sys
from datetime import UTC, datetime
from importlib.metadata import metadata
from pathlib import Path


HERE = Path(__file__).parent


# allow to import our extensions
sys.path[:0] = [str(HERE / 'ext')]


# -- General configuration ------------------------------------------------


# General information
project = 'anndata2ri'
meta = metadata(project)
author = meta['author-email'].split('"')[1]
copyright = f'{datetime.now(tz=UTC):%Y}, {author}.'  # noqa: A001
version = meta['version']
release = version

# default settings
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'

extensions = [
    'scverse_doc',
    'sphinx_autodoc_typehints',
    *[p.stem for p in (HERE / 'ext').glob('*.py')],
]

always_use_bars_union = True
typehints_document_rtype_none = False
typehints_use_rtype = False
autodoc_member_order = 'bysource'
napoleon_custom_sections = ['Supports']

intersphinx_mapping = dict(
    anndata=('https://anndata.scverse.org/en/stable/', None),
    ipython=('https://ipython.readthedocs.io/en/stable/', None),
    numpy=('https://numpy.org/doc/stable/', None),
    pandas=('http://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    rpy2=('https://rpy2.github.io/doc/latest/html/', None),
    scipy=('https://docs.scipy.org/doc/scipy/', None),
)

xref_overrides = {
    # rpy2 documents these modules only as page sections
    ('py:mod', 'rpy2.robjects.numpy2ri'): ('std:ref', 'robjects-numpy'),
    ('py:mod', 'rpy2.robjects.pandas2ri'): ('std:ref', 'robjects-pandas'),
    'rpy2.robjects.conversion.Converter': ('py:mod', 'rpy2.robjects.conversion'),
}


# -- Options for HTML output ----------------------------------------------

source_repository = 'https://github.com/theislab/anndata2ri'
html_static_path = ['_static']
html_css_files = ['css/custom.css']


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f'{project}doc'
doc_title = f'{project} Documentation'
latex_documents = [(master_doc, f'{project}.tex', doc_title, author, 'manual')]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (master_doc, project, doc_title, author, project, 'One line description of project.', 'Miscellaneous')
]
