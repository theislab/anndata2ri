"""Pytest test fixtures and other utilities/enumerations of functionality."""

from __future__ import annotations

from abc import ABC, abstractmethod
from types import ModuleType
from typing import TYPE_CHECKING, Any

import pytest
from rpy2.robjects import default_converter, globalenv
from rpy2.robjects.conversion import Converter, localconverter


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp

    Py2R = Callable[['ConversionModule', Any], Any]
    R2Py = Callable[['ConversionModule', Callable[[], Sexp]], Any]


class ConversionModule(ModuleType, ABC):
    """Type for conversion modules."""

    @property
    @abstractmethod
    def converter(self) -> Converter:
        """Conversion modules have a “converter”."""

    @abstractmethod
    def activate(self) -> None:
        """Conversion modules can be a “activate”d."""

    @abstractmethod
    def deactivate(self) -> None:
        """Conversion modules can be a “deactivate”d."""


def _conversion_py2rpy_manual(conv_mod: ConversionModule, dataset: Any) -> Sexp:  # noqa: ANN401
    converter = conv_mod.get_conversion() if hasattr(conv_mod, 'get_conversion') else conv_mod.converter
    return converter.py2rpy(dataset)


def _conversion_py2rpy_local(conv_mod: ConversionModule, dataset: Any) -> Sexp:  # noqa: ANN401
    converter = conv_mod.get_conversion() if hasattr(conv_mod, 'get_conversion') else conv_mod.converter
    with localconverter(converter):
        globalenv['temp'] = dataset
    return globalenv['temp']


def _conversion_py2rpy_activate(conv_mod: ConversionModule, dataset: Any) -> Sexp:  # noqa: ANN401
    try:
        with pytest.warns(DeprecationWarning, match=r'global conversion'):
            conv_mod.activate()
        globalenv['temp'] = dataset
    finally:
        conv_mod.deactivate()
    return globalenv['temp']


conversions_py2rpy: list[Py2R] = [
    pytest.param(_conversion_py2rpy_manual, id='manual'),
    pytest.param(_conversion_py2rpy_local, id='local'),
    pytest.param(_conversion_py2rpy_activate, id='activate'),
]


@pytest.fixture(params=conversions_py2rpy)
def py2r(request: pytest.FixtureRequest) -> Py2R:
    """Ways to convert a Python object to an R object."""
    return request.param


def _conversion_rpy2py_manual(conv_mod: ConversionModule, dataset: Callable[[], Sexp]) -> Any:  # noqa: ANN401
    converter = conv_mod.get_conversion() if hasattr(conv_mod, 'get_conversion') else conv_mod.converter
    return converter.rpy2py(dataset())


def _conversion_rpy2py_local(conv_mod: ConversionModule, dataset: Callable[[], Sexp]) -> Any:  # noqa: ANN401
    # Needs default_converter to e.g. call `as` on a SummarizedExperiment:
    # Calling a R function returning a S4 object requires py2rpy[RS4], py2rpy[str], …
    converter = conv_mod.get_conversion() if hasattr(conv_mod, 'get_conversion') else conv_mod.converter
    with localconverter(default_converter + converter):
        return dataset()


def _conversion_rpy2py_activate(conv_mod: ConversionModule, dataset: Callable[[], Sexp]) -> Any:  # noqa: ANN401
    try:
        with pytest.warns(DeprecationWarning, match=r'global conversion'):
            conv_mod.activate()
        return dataset()
    finally:
        conv_mod.deactivate()


conversions_rpy2py: list[R2Py] = [
    pytest.param(_conversion_rpy2py_manual, id='manual'),
    pytest.param(_conversion_rpy2py_local, id='local'),
    pytest.param(_conversion_rpy2py_activate, id='activate'),
]


@pytest.fixture(params=conversions_rpy2py)
def r2py(request: pytest.FixtureRequest) -> R2Py:
    """Ways to convert an R object to a Python object."""
    return request.param
