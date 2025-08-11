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

from typing import TYPE_CHECKING

from . import _py2r, _r2py  # noqa: F401
from ._conv import converter as _c
from ._ipython import set_ipython_converter


if TYPE_CHECKING:
    from anndata import AnnData
    from pandas import DataFrame
    from rpy2.rinterface import Sexp


__all__ = ['converter', 'py2rpy', 'rpy2py', 'set_ipython_converter']


converter = _c
"""A converter able to convert most things into an :class:`~anndata.AnnData` object."""


def py2rpy(obj: AnnData) -> Sexp:
    """Convert Python objects to R interface objects.

    Supports
    --------
    :class:`~anndata.AnnData`
        → :rcls:`~SingleCellExperiment::SingleCellExperiment`

    """
    return converter.py2rpy(obj)


def rpy2py(obj: Sexp) -> AnnData | DataFrame:
    """Convert R interface objects to Python objects.

    Supports
    --------
    :rcls:`~SingleCellExperiment::SingleCellExperiment`
        → :class:`~anndata.AnnData`
    :rcls:`S4Vectors::DataFrame`
        → :class:`pandas.DataFrame`

    """
    return converter.rpy2py(obj)
