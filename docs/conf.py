"""Sphinx configuration."""

import sys
from datetime import datetime, timezone
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
copyright = f'{datetime.now(tz=timezone.utc):%Y}, {author}.'  # noqa: A001
version = meta['version']
release = version

# default settings
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary',
    'sphinx_autodoc_typehints',
    'scanpydoc',
    *[p.stem for p in (HERE / 'ext').glob('*.py')],
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
todo_include_todos = False

intersphinx_mapping = dict(
    anndata=('https://anndata.readthedocs.io/en/latest/', None),
    numpy=('https://docs.scipy.org/doc/numpy/', None),
    pandas=('http://pandas.pydata.org/pandas-docs/stable/', None),
    python=('https://docs.python.org/3', None),
    rpy2=('https://rpy2.github.io/doc/latest/html/', None),
    scipy=('https://docs.scipy.org/doc/scipy/reference/', None),
)


# -- Options for HTML output ----------------------------------------------


html_theme = 'scanpydoc'
html_theme_options = dict(collapse_navigation=True)
html_static_path = ['_static']
html_css_files = ['css/custom.css']
html_context = dict(
    display_github=True,
    github_user='theislab',
    github_repo='anndata2ri',
    github_version='main',
    conf_py_path='/docs/',
)


# -- Options for other output formats ------------------------------------------


htmlhelp_basename = f'{project}doc'
doc_title = f'{project} Documentation'
latex_documents = [(master_doc, f'{project}.tex', doc_title, author, 'manual')]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (master_doc, project, doc_title, author, project, 'One line description of project.', 'Miscellaneous')
]
