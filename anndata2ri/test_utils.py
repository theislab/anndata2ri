from abc import ABC, abstractmethod
from types import ModuleType
from typing import List, Callable, Any

import pytest
from rpy2.rinterface import Sexp
from rpy2.robjects import globalenv
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


def conversion_py2rpy_localconv(conv_mod: ConversionModule, dataset: Any) -> Sexp:
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
    pytest.param(conversion_py2rpy_localconv, id="local"),
    pytest.param(conversion_py2rpy_activate, id="activate"),
]
