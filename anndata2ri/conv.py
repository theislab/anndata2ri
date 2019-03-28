from typing import Optional

from rpy2.robjects import conversion, pandas2ri
from rpy2.robjects.conversion import overlay_converter


original_converter: Optional[conversion.Converter] = None
converter = conversion.Converter("original anndata conversion")


def full_converter() -> conversion.Converter:
    pandas2ri.activate()
    new_converter = conversion.Converter("anndata conversion", template=conversion.converter)
    pandas2ri.deactivate()

    overlay_converter(converter, new_converter)

    return new_converter


def activate():
    """
    Activate conversion for :class:`~anndata.AnnData` objects
    as well as :doc:`numpy` arrays and :class:`pandas.DataFrame`\ s
    via ``rpy2.robjects.numpy2ri`` and ``rpy2.robjects.pandas2ri``.
    """
    global original_converter

    if original_converter is not None:
        return

    new_converter = full_converter()
    original_converter = conversion.converter
    conversion.set_conversion(new_converter)


def deactivate():
    """Deactivate!"""
    global original_converter

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
