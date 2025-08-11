from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.core.magic import MagicsManager
    from rpy2.ipython.rmagic import RMagics
    from rpy2.robjects.conversion import Converter


def set_ipython_converter(ipython: InteractiveShell | None = None, converter: Converter | None = None) -> None:
    """Set the default converter for :mod:`~rpy2.ipython.rmagic` in IPython.

    Parameters
    ----------
    ipython
        The IPython instance to set the converter for.
        If not specified, the current IPython instance is used.
    converter
        The converter to use.
        If not specified, :attr:`~anndata2ri.converter` is used.
    """
    import anndata2ri  # noqa: PLC0415

    if ipython is None:
        import IPython  # noqa: PLC0415

        ipython: InteractiveShell | None = IPython.get_ipython()

    if converter is None:
        converter = anndata2ri.converter

    mm: MagicsManager = ipython.magics_manager
    r_magics: RMagics = mm.registry['RMagics']
    r_magics.options.converter = converter
