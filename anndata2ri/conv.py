from typing import Optional

from rpy2.robjects import conversion, numpy2ri, pandas2ri, default_converter, SexpVector
from rpy2.robjects.conversion import overlay_converter

from . import scipy2ri


original_converter: Optional[conversion.Converter] = None
converter = conversion.Converter("original anndata conversion")

mat_converter = default_converter + numpy2ri.converter + scipy2ri.converter
# default_converter has SexpVector registered, so we need to overwrite it.
mat_converter.rpy2py.register(SexpVector, numpy2ri.rpy2py_sexp)


def full_converter() -> conversion.Converter:
    pandas2ri.activate()
    new_converter = conversion.Converter("anndata conversion", template=conversion.converter)
    pandas2ri.deactivate()

    overlay_converter(scipy2ri.converter, new_converter)
    # overwrite the scipy2ri Sexp4 converter and add our others
    overlay_converter(converter, new_converter)

    return new_converter


def activate():
    r"""
    Activate conversion for :class:`~anndata.AnnData` objects
    as well as :doc:`numpy` arrays and :class:`pandas.DataFrame`\ s
    via ``rpy2.robjects.numpy2ri`` and ``rpy2.robjects.pandas2ri``.

    Does nothing if this is the active converter.
    """
    global original_converter

    if original_converter is not None:
        return

    new_converter = full_converter()
    original_converter = conversion.converter
    conversion.set_conversion(new_converter)


def deactivate():
    """Deactivate the conversion described above if it is active."""
    global original_converter

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
