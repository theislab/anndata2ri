import numpy as np
import pytest
import scanpy as sc
from anndata import AnnData
from rpy2.robjects import baseenv, globalenv, r
from rpy2.robjects.conversion import ConversionContext

import anndata2ri


ad_empty = lambda x: None, (0, 0), AnnData
ad_simple = lambda x: None, (2, 3), lambda: AnnData(np.array([[1, 2, 3], [0.3, 0.2, 0.1]]), dict(cluster=[1, 2]))
ad_krumsi = lambda x: None, (640, 11), sc.datasets.krumsiek11


@pytest.mark.parametrize("check,shape,dataset", [ad_empty, ad_simple, ad_krumsi])
def test_py2rpy_manual(check, shape, dataset):
    ex = anndata2ri.converter.py2rpy(dataset())
    assert tuple(baseenv["dim"](ex)[::-1]) == shape
    check(ex)


@pytest.mark.parametrize("check,shape,dataset", [ad_empty, ad_simple, ad_krumsi])
def test_py2rpy_with(check, shape, dataset):
    with ConversionContext(anndata2ri.converter):
        globalenv["adata"] = dataset()
    ex = globalenv["adata"]
    assert tuple(baseenv["dim"](ex)[::-1]) == shape
    check(ex)


@pytest.mark.parametrize("check,shape,dataset", [ad_empty, ad_simple, ad_krumsi])
def test_py2rpy_activate(check, shape, dataset):
    try:
        anndata2ri.activate()
        globalenv["adata"] = dataset()
    finally:
        anndata2ri.deactivate()
    ex = globalenv["adata"]
    assert tuple(baseenv["dim"](ex)[::-1]) == shape
    check(ex)
