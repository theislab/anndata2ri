from typing import Optional

from rpy2.robjects import conversion, numpy2ri

original_converter: Optional[conversion.Converter] = None
converter = conversion.Converter("original scipy conversion")


def activate():
    """Activate!"""
    global original_converter

    if original_converter is not None:
        return

    original_converter = conversion.converter

    numpy2ri.activate()
    new_converter = conversion.Converter("scipy conversion", template=conversion.converter)
    numpy2ri.deactivate()

    conversion.set_conversion(new_converter)


def deactivate():
    """Deactivate!"""
    global original_converter

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
