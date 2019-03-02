import numpy as np
import pytest
import scanpy as sc
from anndata import AnnData
from rpy2.robjects import baseenv, globalenv, r
from rpy2.robjects.conversion import ConversionContext
from rpy2.robjects.packages import importr

import anndata2ri


def mk_ad_simple():
    return AnnData(
        np.array([[1, 2, 3], [0.3, 0.2, 0.1]]),
        dict(cluster=[1, 2]),
        obsm=dict(X_pca=np.array([[1, 3, 5, 7], [2, 4, 6, 8]])),
    )


def check_simple(ex):
    sce = importr("SingleCellExperiment")
    assert [str(n) for n in sce.reducedDimNames(ex)] == ["X_pca"]
    pca = sce.reducedDim(ex, "X_pca")
    assert tuple(baseenv["dim"](pca)) == (2, 4)


ad_empty = lambda x: None, (0, 0), AnnData
ad_simple = check_simple, (2, 3), mk_ad_simple
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
