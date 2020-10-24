r"""
Converter between Python’s AnnData and R’s SingleCellExperiment.


==========================================================  =====  ================================================
:rcls:`~SingleCellExperiment::SingleCellExperiment`                :class:`~anndata.AnnData`
==========================================================  =====  ================================================
:rman:`~SummarizedExperiment::assay`\ ``(d, 'X')``          ⇄      :attr:`~anndata.AnnData.X`
:rman:`~SummarizedExperiment::assay`\ ``(d, 'counts')``     ⇄      :attr:`~anndata.AnnData.layers`\ ``['counts']``
:rman:`~SummarizedExperiment::colData`\ ``(d)``             ⇄      :attr:`~anndata.AnnData.obs`
:rman:`~SummarizedExperiment::rowData`\ ``(d)``             ⇄      :attr:`~anndata.AnnData.var`
:rman:`~S4Vectors::metadata`\ ``(d)``                       ⇄      :attr:`~anndata.AnnData.uns`
:rman:`~SingleCellExperiment::reducedDim`\ ``(d, 'PCA')``   ⇄      :attr:`~anndata.AnnData.obsm`\ ``['X_pca']``
:rman:`~SingleCellExperiment::reducedDim`\ ``(d, 'DM')``    ⇄      :attr:`~anndata.AnnData.obsm`\ ``['X_diffmap']``
==========================================================  =====  ================================================
"""
__all__ = ["activate", "deactivate", "py2rpy", "rpy2py", "converter"]

import traceback
from typing import Any

from get_version import get_version


__author__ = "Philipp Angerer"
__version__ = get_version(__file__)


def within_flit():
    for frame in traceback.extract_stack():
        if frame.name == "get_docstring_and_version_via_import":
            return True
    return False


if not within_flit():
    from rpy2.rinterface import Sexp

    from .conv import converter, activate, deactivate
    from . import py2r, r2py

    def py2rpy(obj: Any) -> Sexp:
        """
        Convert Python objects to R interface objects. Supports:

        - :class:`~anndata.AnnData` → :rcls:`~SingleCellExperiment::SingleCellExperiment`
        """
        return converter.py2rpy(obj)

    def rpy2py(obj: Any) -> Sexp:
        """
        Convert R interface objects to Python objects. Supports:

        - :rcls:`~SingleCellExperiment::SingleCellExperiment` → :class:`~anndata.AnnData`
        - :rcls:`S4Vectors::DataFrame` → :class:`pandas.DataFrame`
        """
        return converter.rpy2py(obj)
