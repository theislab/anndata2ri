import numpy as np
import pytest
from anndata import AnnData
from rpy2.robjects import conversion, baseenv

import anndata2ri


@pytest.fixture
def adata_simple():
    return AnnData(np.array([[1, 2, 3], [0.3, 0.2, 0.1]]), dict(cluster=[1, 2]))


def test_py2rpy_with(adata_simple):
    with conversion.ConversionContext(anndata2ri.create_converter()) as c:
        ex = c.py2rpy(adata_simple)
        assert baseenv["dim"](ex).tolist() == [3, 2]


def test_py2rpy_simple(adata_simple):
    try:
        anndata2ri.activate()
        ex = conversion.py2rpy(adata_simple)
        assert baseenv["dim"](ex).tolist() == [3, 2]
    finally:
        anndata2ri.deactivate()
