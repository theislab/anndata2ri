r"""Convert :mod:`scipy.sparse` matrices between Python and R.

For a detailed comparison between the two languages’
sparse matrix environment, see `issue #8`_.

.. _issue #8: https://github.com/flying-sheep/anndata2ri/issues/8

Here’s an overview over the matching classes (note that ``dtype=float32`` is also supported):

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

from . import _py2r, _r2py  # noqa: F401
from ._conv import converter
from ._support import supported_r_matrix_classes, supported_r_matrix_storage, supported_r_matrix_types


__all__ = [
    'converter',
    'supported_r_matrix_classes',
    'supported_r_matrix_storage',
    'supported_r_matrix_types',
]

converter = converter  # noqa: PLW0127
"""The :class:`~rpy2.robjects.conversion.Converter` for :mod:`scipy.sparse`.

Includes :mod:`rpy2.robjects.numpy2ri`.
"""


supported_r_matrix_types = supported_r_matrix_types  # noqa: PLW0127
"""The Matrix data types supported by this module; Double, Logical, and patterN."""

supported_r_matrix_storage = supported_r_matrix_storage  # noqa: PLW0127
"""The Matrix storage types supported by this module; Column-sparse, Row-Sparse, Triplets, and DIagonal."""
