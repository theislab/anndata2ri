from functools import wraps
from typing import Optional, Any, Callable, Tuple, Type

import numpy as np
from rpy2.rinterface import Sexp
from rpy2.robjects import default_converter, numpy2ri, baseenv
from rpy2.robjects import Vector, BoolVector, IntVector, FloatVector
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import Package
from scipy import sparse

from ..rpy2_ext import importr
from .conv import converter


methods: Optional[Package] = None
as_logical: Optional[Callable[[Any], BoolVector]] = None
as_integer: Optional[Callable[[Any], IntVector]] = None
as_double: Optional[Callable[[Any], FloatVector]] = None


def get_type_conv(dtype: np.dtype) -> Tuple[str, Callable[[np.ndarray], Sexp], Type[Vector]]:
    if np.issubdtype(dtype, np.floating):
        return "d", as_double, FloatVector
    elif np.issubdtype(dtype, np.bool_):
        return "l", as_logical, BoolVector
    else:
        raise ValueError(f"Unknown dtype {dtype!r} cannot be converted to ?gRMatrix.")


def py2r_context(f):
    @wraps(f)
    def wrapper(obj):
        global methods, as_logical, as_integer, as_double
        if methods is None:
            importr("Matrix")  # make class available
            methods = importr("methods")
            as_logical = baseenv["as.logical"]
            as_integer = baseenv["as.integer"]
            as_double = baseenv["as.double"]

        with localconverter(default_converter + numpy2ri.converter):
            return f(obj)

    return wrapper


@converter.py2rpy.register(sparse.csc_matrix)
@py2r_context
def csc_to_rmat(csc: sparse.csc_matrix):
    csc.sort_indices()
    t, conv_data, _ = get_type_conv(csc.dtype)
    return methods.new(
        f"{t}gCMatrix",
        i=as_integer(csc.indices),
        p=as_integer(csc.indptr),
        x=conv_data(csc.data),
        Dim=as_integer(list(csc.shape)),
    )


@converter.py2rpy.register(sparse.csr_matrix)
@py2r_context
def csr_to_rmat(csr: sparse.csr_matrix):
    csr.sort_indices()
    t, conv_data, _ = get_type_conv(csr.dtype)
    return methods.new(
        f"{t}gRMatrix",
        j=as_integer(csr.indices),
        p=as_integer(csr.indptr),
        x=conv_data(csr.data),
        Dim=as_integer(list(csr.shape)),
    )


@converter.py2rpy.register(sparse.coo_matrix)
@py2r_context
def coo_to_rmat(coo: sparse.coo_matrix):
    t, conv_data, _ = get_type_conv(coo.dtype)
    return methods.new(
        f"{t}gTMatrix",
        i=as_integer(coo.row),
        j=as_integer(coo.col),
        x=conv_data(coo.data),
        Dim=as_integer(list(coo.shape)),
    )


@converter.py2rpy.register(sparse.dia_matrix)
@py2r_context
def dia_to_rmat(dia: sparse.dia_matrix):
    t, conv_data, vec_cls = get_type_conv(dia.dtype)
    if len(dia.offsets) > 1:
        raise ValueError(
            "Cannot convert a dia_matrix with more than 1 diagonal to a *diMatrix. "
            f"R diagonal matrices only support 1 diagonal, but this has {len(dia.offsets)}."
        )
    is_unit = np.all(dia.data == 1)
    return methods.new(
        f"{t}diMatrix",
        x=vec_cls([]) if is_unit else conv_data(dia.data),
        diag="U" if is_unit else "N",
        Dim=as_integer(list(dia.shape)),
    )
