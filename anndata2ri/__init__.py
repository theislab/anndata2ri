"""
Converter between Python’s AnnData and R’s SingleCellExperiment.
"""


from get_version import get_version

__version__ = get_version(__file__)


try:  # This is so that flit can import this. There must be a better way.
    from .conv import converter, activate, deactivate
    from . import py2r, r2py

    py2rpy = converter.py2rpy
    rpy2py = converter.rpy2py
except ImportError as e:
    import warnings

    warnings.warn(str(e), ImportWarning)
