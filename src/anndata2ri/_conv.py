from __future__ import annotations

from functools import cache
from typing import TYPE_CHECKING
from uuid import uuid4

import h5py
import numpy as np
import pandas as pd
from anndata.io import sparse_dataset
from rpy2 import rinterface
from rpy2.robjects import conversion, numpy2ri, pandas2ri
from rpy2.robjects.vectors import IntVector, StrVector

from . import scipy2ri


if TYPE_CHECKING:
    from anndata.abc import CSCDataset, CSRDataset
    from rpy2.rinterface import Sexp
    from scipy.sparse import sparray, spmatrix


converter = conversion.Converter(
    'original anndata conversion', template=numpy2ri.converter + pandas2ri.converter + scipy2ri.converter
)


@converter.py2rpy.register(np.ndarray)
def py2rpy_ndarray(obj: np.ndarray) -> Sexp:
    # numpy2ri doesn’t know numpy’s variable-width StringDType (kind 'T') yet
    if obj.dtype.kind != 'T':
        return numpy2ri.numpy2rpy(obj)
    obj = obj.astype(np.dtypes.StringDType(na_object=rinterface.NA_Character))
    res = StrVector(obj.ravel(order='F').tolist())
    res.do_slot_assign('dim', IntVector(obj.shape))  # like numpy2ri does
    return res


_mat_converter = numpy2ri.converter + scipy2ri.converter


def mat_py2rpy(obj: np.ndarray | sparray | spmatrix | pd.DataFrame) -> Sexp:
    if isinstance(obj, pd.DataFrame):
        numeric_cols = obj.dtypes.map(pd.api.types.is_numeric_dtype)
        if not numeric_cols.all():
            non_num = numeric_cols.index[~numeric_cols]
            msg = f'DataFrame contains non-numeric columns {list(non_num)}'
            raise ValueError(msg)
        obj = obj.to_numpy()
    return _mat_converter.py2rpy(obj)


@cache
def _r_realize_delayed() -> rinterface.SexpClosure:
    return rinterface.evalr("""function(x) {
        if (!is(x, "DelayedArray")) NULL
        else if (DelayedArray::is_sparse(x)) as(x, "CsparseMatrix")
        else as.matrix(x)
    }""")


def realize_delayed(obj: Sexp) -> Sexp | None:
    """Load a DelayedArray (e.g. BiocSingular’s LowRankMatrix) into memory, or return `None` if it isn’t one."""
    if not isinstance(obj, rinterface.SexpS4):
        return None
    mat = _r_realize_delayed()(obj)
    return None if isinstance(mat, rinterface.NULLType) else mat


@cache
def _r_h5_location() -> rinterface.SexpClosure:
    # A CSC matrix in R (genes × cells) is a CSR matrix from Python’s (cells × genes) point of view and vice versa
    return rinterface.evalr("""function(x) {
        if (!is(x, "DelayedArray") || !DelayedArray::isPristine(x)) return(NULL)
        s <- DelayedArray::seed(x)
        sparse <- function(format) list(
            format = format, path = s@filepath, name = s@group,
            data = if (is.null(s@subdata)) "data" else s@subdata, shape = rev(dim(s))
        )
        if (is(s, "HDF5ArraySeed")) list(format = "dense", path = s@filepath, name = s@name)
        else if (is(s, "CSC_H5SparseMatrixSeed")) sparse("csr")
        else if (is(s, "CSR_H5SparseMatrixSeed")) sparse("csc")
        else NULL
    }""")


def _open_h5(loc: rinterface.ListSexpVector) -> h5py.Dataset | CSRDataset | CSCDataset:
    info = dict(zip(loc.names, loc, strict=True))
    fmt, path, name = info['format'][0], info['path'][0], info['name'][0]
    if fmt == 'dense':
        return h5py.File(path, 'r')[name]
    # R’s sparse groups lack anndata’s encoding metadata, so link them into an in-memory file that has it
    # HDF5 tracks open files by name even when they’re in memory, so make it unique
    wrapper = h5py.File(f'{path}/{name}/{uuid4()}', 'w', driver='core', backing_store=False)
    for key, src in [('data', info['data'][0]), ('indices', 'indices'), ('indptr', 'indptr')]:
        wrapper[key] = h5py.ExternalLink(path, f'{name}/{src}')
    wrapper.attrs.update({'encoding-type': f'{fmt}_matrix', 'encoding-version': '0.1.0', 'shape': tuple(info['shape'])})
    return sparse_dataset(wrapper)


def mat_rpy2py(
    obj: Sexp, *, transpose: bool = False
) -> np.ndarray | sparray | h5py.Dataset | CSRDataset | CSCDataset | Sexp:
    """Convert an R matrix to Python, optionally transposing it.

    When transposing, pristine HDF5-backed DelayedArrays are returned as backed arrays,
    since their on-disk layout is already transposed from R’s point of view.
    """
    if transpose and isinstance(obj, rinterface.SexpS4):
        loc = _r_h5_location()(obj)
        if not isinstance(loc, rinterface.NULLType):
            return _open_h5(loc)
    if (mat := realize_delayed(obj)) is not None:
        obj = mat
    conv = _mat_converter.rpy2py(obj)
    return conv.T if transpose and not isinstance(conv, rinterface.Sexp) else conv
