from __future__ import annotations

from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
import pandas as pd
from rpy2.robjects import conversion, numpy2ri, pandas2ri

from . import scipy2ri


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp
    from scipy.sparse import spmatrix


original_converter: conversion.Converter | None = None
converter = (
    conversion.Converter('original anndata conversion', template=conversion.get_conversion())
    + numpy2ri.converter
    + pandas2ri.converter
    + scipy2ri.converter
)


_mat_converter = numpy2ri.converter + scipy2ri.converter


def mat_py2rpy(obj: np.ndarray | spmatrix | pd.DataFrame) -> Sexp:
    if isinstance(obj, pd.DataFrame):
        numeric_cols = obj.dtypes.map(lambda dt: np.issubdtype(dt, np.number))
        if not numeric_cols.all():
            non_num = numeric_cols.index[~numeric_cols]
            msg = f'DataFrame contains non-numeric columns {list(non_num)}'
            raise ValueError(msg)
        obj = obj.to_numpy()
    return _mat_converter.py2rpy(obj)


mat_rpy2py: Callable[[Sexp], np.ndarray | spmatrix | Sexp] = _mat_converter.rpy2py


def activate() -> None:
    r"""Activate conversion for supported objects.

    This includes :class:`~anndata.AnnData` objects
    as well as :ref:`numpy:arrays` and :class:`pandas.DataFrame`\ s
    via ``rpy2.robjects.numpy2ri`` and ``rpy2.robjects.pandas2ri``.

    Does nothing if this is the active converter.
    """
    warn(
        'The global conversion available with activate() '
        'is deprecated and will be removed in the next major release. '
        'Use a local converter.',
        category=DeprecationWarning,
        stacklevel=2,
    )
    global original_converter  # noqa: PLW0603

    if original_converter is not None:
        return

    new_converter = converter
    original_converter = conversion.get_conversion()
    conversion.set_conversion(new_converter)


def deactivate() -> None:
    """Deactivate the conversion described above if it is active."""
    global original_converter  # noqa: PLW0603

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
