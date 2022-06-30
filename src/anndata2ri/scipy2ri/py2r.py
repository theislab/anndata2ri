from functools import wraps
from typing import Callable, Optional

import numpy as np
from rpy2.rinterface import Sexp
from rpy2.robjects import default_converter, numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import Package, SignatureTranslatedAnonymousPackage
from scipy import sparse

from ..rpy2_ext import importr
from .conv import converter


matrix: Optional[SignatureTranslatedAnonymousPackage] = None
base: Optional[Package] = None


def get_type_conv(dtype: np.dtype) -> Callable[[np.ndarray], Sexp]:
    global base
    if base is None:
        base = importr('base')
    if np.issubdtype(dtype, np.floating):
        return base.as_double
    elif np.issubdtype(dtype, np.bool_):
        return base.as_logical
    else:
        raise ValueError(f'Unknown dtype {dtype!r} cannot be converted to ?gRMatrix.')


def py2r_context(f):
    @wraps(f)
    def wrapper(obj):
        global as_logical, as_integer, as_double, matrix
        if matrix is None:
            importr('Matrix')  # make class available
            matrix = SignatureTranslatedAnonymousPackage(
                """
                sparse_matrix <- function(x, conv_data, dims, ...) {
                    Matrix::sparseMatrix(
                        ...,
                        x=conv_data(x),
                        dims=as.integer(dims),
                        index1=FALSE
                    )
                }

                from_csc <- function(i, p, x, dims, conv_data) {
                    sparse_matrix(
                        i=as.integer(i),
                        p=as.integer(p),
                        x=x,
                        conv_data=conv_data,
                        dims=dims,
                        repr="C"
                    )
                }

                from_csr <- function(j, p, x, dims, conv_data) {
                    sparse_matrix(
                        j=as.integer(j),
                        p=as.integer(p),
                        x=x,
                        conv_data=conv_data,
                        dims=dims,
                        repr="R"
                    )
                }

                from_coo <- function(i, j, x, dims, conv_data) {
                    sparse_matrix(
                        i=as.integer(i),
                        j=as.integer(j),
                        x=x,
                        conv_data=conv_data,
                        dims=dims,
                        repr="T"
                    )
                }

                from_dia <- function(n, x, conv_data) {
                    Matrix::Diagonal(n=as.integer(n), x=conv_data(x))
                }
                """,
                'matrix',
            )

        return f(obj)

    return wrapper


@converter.py2rpy.register(sparse.csc_matrix)
@py2r_context
def csc_to_rmat(csc: sparse.csc_matrix):
    csc.sort_indices()
    conv_data = get_type_conv(csc.dtype)
    with localconverter(default_converter + numpy2ri.converter):
        return matrix.from_csc(i=csc.indices, p=csc.indptr, x=csc.data, dims=list(csc.shape), conv_data=conv_data)


@converter.py2rpy.register(sparse.csr_matrix)
@py2r_context
def csr_to_rmat(csr: sparse.csr_matrix):
    csr.sort_indices()
    conv_data = get_type_conv(csr.dtype)
    with localconverter(default_converter + numpy2ri.converter):
        return matrix.from_csr(
            j=csr.indices,
            p=csr.indptr,
            x=csr.data,
            conv_data=conv_data,
            dims=list(csr.shape),
        )


@converter.py2rpy.register(sparse.coo_matrix)
@py2r_context
def coo_to_rmat(coo: sparse.coo_matrix):
    conv_data = get_type_conv(coo.dtype)
    with localconverter(default_converter + numpy2ri.converter):
        return matrix.from_coo(
            i=coo.row,
            j=coo.col,
            x=coo.data,
            conv_data=conv_data,
            dims=list(coo.shape),
        )


@converter.py2rpy.register(sparse.dia_matrix)
@py2r_context
def dia_to_rmat(dia: sparse.dia_matrix):
    conv_data = get_type_conv(dia.dtype)
    if len(dia.offsets) > 1:
        raise ValueError(
            'Cannot convert a dia_matrix with more than 1 diagonal to a *diMatrix. '
            f'R diagonal matrices only support 1 diagonal, but this has {len(dia.offsets)}.'
        )
    with localconverter(default_converter + numpy2ri.converter):
        return matrix.from_dia(
            n=dia.shape[0],
            x=dia.data,
            conv_data=conv_data,
        )
