from __future__ import annotations

import numpy as np
import pandas as pd
from rpy2.robjects import conversion, numpy2ri, pandas2ri
from rpy2.robjects.conversion import overlay_converter

from . import scipy2ri


original_converter: conversion.Converter | None = None
converter = conversion.Converter('original anndata conversion')

_mat_converter = numpy2ri.converter + scipy2ri.converter


def mat_py2rpy(obj: np.ndarray) -> np.ndarray:
    if isinstance(obj, pd.DataFrame):
        numeric_cols = obj.dtypes <= np.number
        if not numeric_cols.all():
            raise ValueError('DataFrame contains non-numeric columns')
        obj = obj.to_numpy()
    return _mat_converter.py2rpy(obj)


mat_rpy2py = _mat_converter.rpy2py


def full_converter() -> conversion.Converter:
    pandas2ri.activate()
    new_converter = conversion.Converter('anndata conversion', template=conversion.get_conversion())
    pandas2ri.deactivate()

    overlay_converter(scipy2ri.converter, new_converter)
    # overwrite the scipy2ri Sexp4 converter and add our others
    overlay_converter(converter, new_converter)

    return new_converter


def activate():
    r"""
    Activate conversion for :class:`~anndata.AnnData` objects
    as well as :ref:`numpy:arrays` and :class:`pandas.DataFrame`\ s
    via ``rpy2.robjects.numpy2ri`` and ``rpy2.robjects.pandas2ri``.

    Does nothing if this is the active converter.
    """
    global original_converter

    if original_converter is not None:
        return

    new_converter = full_converter()
    original_converter = conversion.get_conversion()
    conversion.set_conversion(new_converter)


def deactivate():
    """Deactivate the conversion described above if it is active."""
    global original_converter

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
