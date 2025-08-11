from __future__ import annotations

from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.core.magic import MagicsManager
    from rpy2.ipython.rmagic import RMagics
    from rpy2.robjects.conversion import Converter


def set_ipython_converter(ipython: InteractiveShell | None = None, converter: Converter | None = None) -> None:
    import anndata2ri  # noqa: PLC0415

    if ipython is None:
        import IPython  # noqa: PLC0415

        ipython: InteractiveShell | None = IPython.get_ipython()

    if converter is None:
        converter = anndata2ri.converter

    mm: MagicsManager = ipython.magics_manager
    r_magics: RMagics = mm.registry['RMagics']
    r_magics.options.converter = converter
