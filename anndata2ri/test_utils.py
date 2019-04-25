from abc import ABC, abstractmethod
from types import ModuleType
from typing import List, Callable, Any

import pytest
from rpy2.rinterface import Sexp
from rpy2.robjects import globalenv, default_converter
from rpy2.robjects.conversion import localconverter, Converter


class ConversionModule(ModuleType, ABC):
    @property
    @abstractmethod
    def converter(self) -> Converter:
        pass

    @abstractmethod
    def activate(self) -> None:
        pass

    @abstractmethod
    def deactivate(self) -> None:
        pass


def conversion_py2rpy_manual(conv_mod: ConversionModule, dataset: Any) -> Sexp:
    return conv_mod.converter.py2rpy(dataset)


def conversion_py2rpy_local(conv_mod: ConversionModule, dataset: Any) -> Sexp:
    with localconverter(conv_mod.converter):
        globalenv["temp"] = dataset
    return globalenv["temp"]


def conversion_py2rpy_activate(conv_mod: ConversionModule, dataset: Any) -> Sexp:
    try:
        conv_mod.activate()
        globalenv["temp"] = dataset
    finally:
        conv_mod.deactivate()
    return globalenv["temp"]


conversions_py2rpy: List[Callable[[ConversionModule, Any], Sexp]] = [
    pytest.param(conversion_py2rpy_manual, id="manual"),
    pytest.param(conversion_py2rpy_local, id="local"),
    pytest.param(conversion_py2rpy_activate, id="activate"),
]


def conversion_rpy2py_manual(conv_mod: ConversionModule, dataset: Callable[[], Sexp]) -> Any:
    return conv_mod.converter.rpy2py(dataset())


def conversion_rpy2py_local(conv_mod: ConversionModule, dataset: Callable[[], Sexp]) -> Any:
    # Needs default_converter to e.g. call `as` on a SummarizedExperiment:
    # Calling a R function returning a S4 object requires py2rpy[RS4], py2rpy[str], …
    with localconverter(default_converter + conv_mod.converter):
        return dataset()


def conversion_rpy2py_activate(conv_mod: ConversionModule, dataset: Callable[[], Sexp]) -> Any:
    try:
        conv_mod.activate()
        return dataset()
    finally:
        conv_mod.deactivate()


conversions_rpy2py: List[Callable[[ConversionModule, Callable[[], Sexp]], Any]] = [
    pytest.param(conversion_rpy2py_manual, id="manual"),
    pytest.param(conversion_rpy2py_local, id="local"),
    pytest.param(conversion_rpy2py_activate, id="activate"),
]
