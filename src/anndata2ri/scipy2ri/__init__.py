r"""Convert scipy.sparse matrices between Python and R.

For a detailed comparison between the two languages’
sparse matrix environment, see `issue #8`_.

.. _issue #8: https://github.com/flying-sheep/anndata2ri/issues/8

Here’s an overview over the matching classes

=====================================================  ======================================================
R                                                      Python
=====================================================  ======================================================
:rcls:`~Matrix::dgCMatrix`                             :class:`~scipy.sparse.csc_matrix`\ ``(dtype=float64)``
:rcls:`~Matrix::lgCMatrix`/:rcls:`~Matrix::ngCMatrix`  :class:`~scipy.sparse.csc_matrix`\ ``(dtype=bool)``
:rcls:`~Matrix::dgRMatrix`                             :class:`~scipy.sparse.csr_matrix`\ ``(dtype=float64)``
:rcls:`~Matrix::lgRMatrix`/:rcls:`~Matrix::ngRMatrix`  :class:`~scipy.sparse.csr_matrix`\ ``(dtype=bool)``
:rcls:`~Matrix::dgTMatrix`                             :class:`~scipy.sparse.coo_matrix`\ ``(dtype=float64)``
:rcls:`~Matrix::lgTMatrix`/:rcls:`~Matrix::ngTMatrix`  :class:`~scipy.sparse.coo_matrix`\ ``(dtype=bool)``
:rcls:`~Matrix::ddiMatrix`                             :class:`~scipy.sparse.dia_matrix`\ ``(dtype=float64)``
:rcls:`~Matrix::ldiMatrix`                             :class:`~scipy.sparse.dia_matrix`\ ``(dtype=bool)``
=====================================================  ======================================================
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from . import _py2r, _r2py  # noqa: F401
from ._conv import activate, converter, deactivate
from ._support import supported_r_matrix_classes, supported_r_matrix_storage, supported_r_matrix_types


if TYPE_CHECKING:
    from rpy2.rinterface import Sexp
    from scipy import sparse


__all__ = [
    'activate',
    'deactivate',
    'py2rpy',
    'rpy2py',
    'converter',
    'supported_r_matrix_types',
    'supported_r_matrix_storage',
    'supported_r_matrix_classes',
]


supported_r_matrix_types = supported_r_matrix_types  # noqa: PLW0127
"""The Matrix data types supported by this module; Double, Logical, and patterN."""

supported_r_matrix_storage = supported_r_matrix_storage  # noqa: PLW0127
"""The Matrix storage types supported by this module; Column-sparse, Row-Sparse, Triplets, and DIagonal."""


def py2rpy(obj: sparse.spmatrix) -> Sexp:
    """Convert scipy sparse matrices objects to R sparse matrices.

    Supports:

    :class:`~scipy.sparse.csc_matrix` (dtype in {float32, float64, bool}) →
        :rcls:`~Matrix::dgCMatrix` or :rcls:`~Matrix::lgCMatrix`
    :class:`~scipy.sparse.csr_matrix` (dtype in {float32, float64, bool}) →
        :rcls:`~Matrix::dgRMatrix` or :rcls:`~Matrix::lgRMatrix`
    :class:`~scipy.sparse.coo_matrix` (dtype in {float32, float64, bool}) →
        :rcls:`~Matrix::dgTMatrix` or :rcls:`~Matrix::lgTMatrix`
    :class:`~scipy.sparse.dia_matrix` (dtype in {float32, float64, bool}) →
        :rcls:`~Matrix::ddiMatrix` or :rcls:`~Matrix::ldiMatrix`
    """
    return converter.py2rpy(obj)


def rpy2py(obj: Sexp) -> sparse.spmatrix:
    """Convert R sparse matrices to scipy sparse matrices.

    Supports:

    :rcls:`~Matrix::dgCMatrix`, :rcls:`~Matrix::lgCMatrix`, or :rcls:`~Matrix::ngCMatrix` →
        :class:`~scipy.sparse.csc_matrix` (dtype float64 or bool)
    :rcls:`~Matrix::dgRMatrix`, :rcls:`~Matrix::lgRMatrix`, or :rcls:`~Matrix::ngRMatrix` →
        :class:`~scipy.sparse.csr_matrix` (dtype float64 or bool)
    :rcls:`~Matrix::dgTMatrix`, :rcls:`~Matrix::lgTMatrix`, or :rcls:`~Matrix::ngTMatrix` →
        :class:`~scipy.sparse.coo_matrix` (dtype float64 or bool)
    :rcls:`~Matrix::ddiMatrix` or :rcls:`~Matrix::ldiMatrix` →
        :class:`~scipy.sparse.dia_matrix` (dtype float64 or bool)
    """
    return converter.rpy2py(obj)
