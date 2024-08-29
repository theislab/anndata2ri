from __future__ import annotations

from functools import cache, wraps
from importlib.resources import files
from typing import TYPE_CHECKING

import numpy as np
from rpy2.robjects import default_converter, numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import InstalledSTPackage, SignatureTranslatedAnonymousPackage
from scipy import sparse

from anndata2ri._rpy2_ext import R_INT_BYTES, importr

from ._conv import converter


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp


@cache
def baseenv() -> InstalledSTPackage:
    return importr('base')


@cache
def matrixenv() -> SignatureTranslatedAnonymousPackage:
    importr('Matrix')  # make class available
    r_code = files('anndata2ri').joinpath('scipy2ri', '_py2r_helpers.r').read_text()
    return SignatureTranslatedAnonymousPackage(r_code, 'matrix')


def get_type_conv(dtype: np.dtype) -> Callable[[np.ndarray], Sexp]:
    base = baseenv()
    if np.issubdtype(dtype, np.floating):
        return base.as_double
    if np.issubdtype(dtype, np.integer):
        if dtype.itemsize <= R_INT_BYTES:
            return base.as_integer
        return base.as_numeric  # maybe uses R_xlen_t?
    if np.issubdtype(dtype, np.bool_):
        return base.as_logical
    msg = f'Unknown dtype {dtype!r} cannot be converted to ?gRMatrix.'
    raise ValueError(msg)


def py2r_context(f: Callable[[sparse.spmatrix], Sexp]) -> Callable[[sparse.spmatrix], Sexp]:
    """R globalenv context with some helper functions."""

    @wraps(f)
    def wrapper(obj: sparse.spmatrix) -> Sexp:
        matrixenv()  # make Matrix class available
        return f(obj)

    return wrapper


@converter.py2rpy.register(sparse.csc_matrix)
@py2r_context
def csc_to_rmat(csc: sparse.csc_matrix) -> Sexp:
    csc.sort_indices()
    conv_data = get_type_conv(csc.dtype)
    with localconverter(default_converter + numpy2ri.converter):
        return matrixenv().from_csc(i=csc.indices, p=csc.indptr, x=csc.data, dims=list(csc.shape), conv_data=conv_data)


@converter.py2rpy.register(sparse.csr_matrix)
@py2r_context
def csr_to_rmat(csr: sparse.csr_matrix) -> Sexp:
    csr.sort_indices()
    conv_data = get_type_conv(csr.dtype)
    with localconverter(default_converter + numpy2ri.converter):
        return matrixenv().from_csr(
            j=csr.indices,
            p=csr.indptr,
            x=csr.data,
            conv_data=conv_data,
            dims=list(csr.shape),
        )


@converter.py2rpy.register(sparse.coo_matrix)
@py2r_context
def coo_to_rmat(coo: sparse.coo_matrix) -> Sexp:
    conv_data = get_type_conv(coo.dtype)
    with localconverter(default_converter + numpy2ri.converter):
        return matrixenv().from_coo(
            i=coo.row,
            j=coo.col,
            x=coo.data,
            conv_data=conv_data,
            dims=list(coo.shape),
        )


@converter.py2rpy.register(sparse.dia_matrix)
@py2r_context
def dia_to_rmat(dia: sparse.dia_matrix) -> Sexp:
    conv_data = get_type_conv(dia.dtype)
    if len(dia.offsets) > 1:
        msg = (
            'Cannot convert a dia_matrix with more than 1 diagonal to a *diMatrix. '
            f'R diagonal matrices only support 1 diagonal, but this has {len(dia.offsets)}.'
        )
        raise ValueError(msg)
    with localconverter(default_converter + numpy2ri.converter):
        return matrixenv().from_dia(
            n=dia.shape[0],
            x=dia.data,
            conv_data=conv_data,
        )
