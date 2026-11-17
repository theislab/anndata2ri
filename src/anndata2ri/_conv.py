from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from rpy2 import rinterface
from rpy2.robjects import conversion, numpy2ri, pandas2ri
from rpy2.robjects.vectors import IntVector, StrVector

from . import scipy2ri


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp
    from scipy.sparse import spmatrix


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


def mat_py2rpy(obj: np.ndarray | spmatrix | pd.DataFrame) -> Sexp:
    if isinstance(obj, pd.DataFrame):
        numeric_cols = obj.dtypes.map(pd.api.types.is_numeric_dtype)
        if not numeric_cols.all():
            non_num = numeric_cols.index[~numeric_cols]
            msg = f'DataFrame contains non-numeric columns {list(non_num)}'
            raise ValueError(msg)
        obj = obj.to_numpy()
    return _mat_converter.py2rpy(obj)


mat_rpy2py: Callable[[Sexp], np.ndarray | spmatrix | Sexp] = _mat_converter.rpy2py
