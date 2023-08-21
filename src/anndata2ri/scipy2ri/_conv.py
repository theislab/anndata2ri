from __future__ import annotations

from warnings import warn

from rpy2.robjects import conversion, numpy2ri
from rpy2.robjects.conversion import overlay_converter


original_converter: conversion.Converter | None = None
converter = conversion.Converter('original scipy conversion')


def activate() -> None:
    """Activate conversion between sparse matrices from Scipy and R’s Matrix package.

    Does nothing if this is the active conversion.
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

    original_converter = conversion.get_conversion()

    numpy2ri.activate()
    new_converter = conversion.Converter('scipy conversion', template=conversion.get_conversion())
    numpy2ri.deactivate()

    overlay_converter(converter, new_converter)

    conversion.set_conversion(new_converter)


def deactivate() -> None:
    """Deactivate the conversion described above if it is active."""
    global original_converter  # noqa: PLW0603

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
