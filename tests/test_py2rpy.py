import numpy as np
import pytest
from anndata import AnnData
from rpy2.robjects import baseenv, globalenv, r
from rpy2.robjects.conversion import ConversionContext

import anndata2ri


@pytest.fixture
def adata_simple():
    return AnnData(np.array([[1, 2, 3], [0.3, 0.2, 0.1]]), dict(cluster=[1, 2]))


def test_py2rpy_manual(adata_simple):
    ex = anndata2ri.converter.py2rpy(adata_simple)
    assert list(baseenv["dim"](ex)) == [3, 2]


def test_py2rpy_with(adata_simple):
    with ConversionContext(anndata2ri.converter):
        globalenv["adata"] = adata_simple
    assert list(r("dim(adata)")) == [3, 2]


def test_py2rpy_activate(adata_simple):
    try:
        anndata2ri.activate()
        globalenv["adata"] = adata_simple
    finally:
        anndata2ri.deactivate()
    assert list(r("dim(adata)")) == [3, 2]
