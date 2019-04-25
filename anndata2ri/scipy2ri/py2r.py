from functools import wraps

import numpy as np
from rpy2.robjects import default_converter, numpy2ri  # , baseenv as base
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from scipy import sparse

from .conv import converter


def get_type_letter(dtype: np.dtype):
    if np.issubdtype(dtype, np.floating):
        return "d"
    elif np.issubdtype(dtype, bool):
        return "l"
    else:
        raise ValueError(f"Unknown dtype {dtype!r} cannot be converted to ?gRMatrix.")


def py2r_context(f):
    @wraps(f)
    def wrapper(obj):
        importr("Matrix")  # make class available
        with localconverter(default_converter + numpy2ri.converter):
            return f(obj)

    return wrapper


@converter.rpy2py.register(sparse.csr_matrix)
@py2r_context
def csr_to_rmat(csr: sparse.csr_matrix):
    methods = importr("methods")
    t = get_type_letter(csr.dtype)
    return methods.new(f"{t}gRMatrix", i=csr.indices, p=csr.indptr, x=csr.data, Dim=list(csr.shape))


@converter.rpy2py.register(sparse.csc_matrix)
@py2r_context
def csc_to_rmat(csc: sparse.csc_matrix):
    methods = importr("methods")
    t = get_type_letter(csc.dtype)
    return methods.new(f"{t}gCMatrix", j=csc.indices, p=csc.indptr, x=csc.data, Dim=list(csc.shape))


@converter.rpy2py.register(sparse.coo_matrix)
@py2r_context
def coo_to_rmat(coo: sparse.coo_matrix):
    methods = importr("methods")
    t = get_type_letter(coo.dtype)
    return methods.new(f"{t}gTMatrix", i=coo.row, j=coo.col, x=coo.data, Dim=list(coo.shape))


@converter.rpy2py.register(sparse.dia_matrix)
@py2r_context
def dia_to_rmat(dia: sparse.dia_matrix):
    methods = importr("methods")
    t = get_type_letter(dia.dtype)
    if len(dia.offsets) > 1:
        raise ValueError(
            "Cannot convert a dia_matrix with more than 1 diagonal to a *diMatrix. "
            f"R diagonal matrices only support 1 diagonal, but this has {len(dia.offsets)}."
        )
    return methods.new(f"{t}gTMatrix", x=dia.data, diag="U" if np.all(dia.data == 1) else "N", Dim=list(dia.shape))
