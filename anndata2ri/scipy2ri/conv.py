from typing import Optional

from rpy2.robjects import conversion, numpy2ri
from rpy2.robjects.conversion import overlay_converter

original_converter: Optional[conversion.Converter] = None
converter = conversion.Converter("original scipy conversion")


def activate():
    """
    Activate conversion between sparse matrices from Scipy and R’s Matrix package.

    Does nothing if this is the active conversion.
    """
    global original_converter

    if original_converter is not None:
        return

    original_converter = conversion.converter

    numpy2ri.activate()
    new_converter = conversion.Converter("scipy conversion", template=conversion.converter)
    numpy2ri.deactivate()

    overlay_converter(converter, new_converter)

    conversion.set_conversion(new_converter)


def deactivate():
    """Deactivate the conversion described above if it is active."""
    global original_converter

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
