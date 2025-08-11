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

from . import _py2r, _r2py  # noqa: F401
from ._conv import converter
from ._ipython import set_ipython_converter


__all__ = ['converter', 'set_ipython_converter']


converter = converter  # noqa: PLW0127
"""A :class:`~rpy2.robjects.conversion.Converter` for :mod:`anndata`.

Includes :mod:`rpy2.robjects.numpy2ri`, :mod:`rpy2.robjects.pandas2ri`, and :mod:`.scipy2ri`.
"""
