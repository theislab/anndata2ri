r"""Converter between Python’s AnnData and R’s SingleCellExperiment.

==========================================================  =  ========================================================
:rcls:`~SingleCellExperiment::SingleCellExperiment`            :class:`~anndata.AnnData`
==========================================================  =  ========================================================
:rman:`~SummarizedExperiment::assay`\ ``(d, 'X')``          ⇄  ``d.``\ :attr:`~anndata.AnnData.X`
:rman:`~SummarizedExperiment::assay`\ ``(d, 'counts')``     ⇄  ``d.``\ :attr:`~anndata.AnnData.layers`\ ``['counts']``
:rman:`~SummarizedExperiment::colData`\ ``(d)``             ⇄  ``d.``\ :attr:`~anndata.AnnData.obs`
:rman:`~SummarizedExperiment::rowData`\ ``(d)``             ⇄  ``d.``\ :attr:`~anndata.AnnData.var`
:rman:`~S4Vectors::metadata`\ ``(d)``                       ⇄  ``d.``\ :attr:`~anndata.AnnData.uns`
:rman:`~SingleCellExperiment::reducedDim`\ ``(d, 'PCA')``   ⇄  ``d.``\ :attr:`~anndata.AnnData.obsm`\ ``['X_pca']``
:rman:`~SingleCellExperiment::reducedDim`\ ``(d, 'DM')``    ⇄  ``d.``\ :attr:`~anndata.AnnData.obsm`\ ``['X_diffmap']``
==========================================================  =  ========================================================
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from . import _py2r, _r2py  # noqa: F401
from ._conv import activate, converter, deactivate


if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame
    from rpy2.rinterface import Sexp


__all__ = ['__version__', 'activate', 'deactivate', 'py2rpy', 'rpy2py', 'converter']

HERE = Path(__file__).parent

__author__ = 'Philipp Angerer'
try:
    from setuptools_scm import get_version

    __version__ = get_version(str(HERE.parent.parent))
except (ImportError, LookupError):
    try:
        from ._version import __version__
    except ImportError as e:
        msg = 'Cannot infer version. Make sure to `pip install` the project or install `setuptools-scm`.'
        raise ImportError(msg) from e


def py2rpy(obj: AnnData) -> Sexp:
    """Convert Python objects to R interface objects.

    Supports:

    - :class:`~anndata.AnnData` → :rcls:`~SingleCellExperiment::SingleCellExperiment`
    """
    return converter.py2rpy(obj)


def rpy2py(obj: Sexp) -> AnnData | DataFrame:
    """Convert R interface objects to Python objects.

    Supports:

    - :rcls:`~SingleCellExperiment::SingleCellExperiment` → :class:`~anndata.AnnData`
    - :rcls:`S4Vectors::DataFrame` → :class:`pandas.DataFrame`
    """
    return converter.rpy2py(obj)
