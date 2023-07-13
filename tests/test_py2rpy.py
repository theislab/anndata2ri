from warnings import catch_warnings, simplefilter

import numpy as np
import pytest
import scanpy as sc
from anndata import AnnData
from pandas import DataFrame
from rpy2.robjects import baseenv, globalenv
from rpy2.robjects.conversion import localconverter

import anndata2ri
from anndata2ri.rpy2_ext import importr


def mk_ad_simple():
    return AnnData(
        np.array([[1, 2, 3], [0.3, 0.2, 0.1]]),
        dict(cluster=[1, 2]),
        obsm=dict(X_pca=np.array([[1, 3, 5, 7], [2, 4, 6, 8]])),
    )


def check_empty(ex):
    pass


def check_pca(ex):
    sce = importr('SingleCellExperiment')
    assert [str(n) for n in sce.reducedDimNames(ex)] == ['PCA']
    pca = sce.reducedDim(ex, 'PCA')
    assert tuple(baseenv['dim'](pca)) == (2, 4)


datasets = [
    pytest.param(check_empty, (0, 0), AnnData, id='empty'),
    pytest.param(check_pca, (2, 3), mk_ad_simple, id='simple'),
    pytest.param(check_empty, (640, 11), sc.datasets.krumsiek11, id='krumsiek'),
    # pytest.param(check_empty, (2730, 3451), sc.datasets.paul15, id="paul"),
]


@pytest.mark.parametrize('check,shape,dataset', datasets)
def test_py2rpy(py2r, check, shape, dataset):
    if dataset is sc.datasets.krumsiek11:
        with pytest.warns(UserWarning, match=r'Duplicated obs_names'):
            ex = py2r(anndata2ri, dataset())
    else:
        ex = py2r(anndata2ri, dataset())
    assert tuple(baseenv['dim'](ex)[::-1]) == shape
    check(ex)


def test_py2rpy2_numpy_pbmc68k():
    """This has some weird metadata"""
    from scanpy.datasets import pbmc68k_reduced

    try:
        anndata2ri.activate()
        with catch_warnings(record=True) as logs:
            simplefilter('ignore', DeprecationWarning)
            globalenv['adata'] = pbmc68k_reduced()
        assert len(logs) == 0, [m.message for m in logs]
    finally:
        anndata2ri.deactivate()


def test_obsm_df():
    """Obsm can contain dataframes"""
    adata = mk_ad_simple()
    adata.obsm['a'] = DataFrame([[5, 6, 7], [8, 9, 0]], index=adata.obs_names)

    with localconverter(anndata2ri.converter):
        globalenv['adata_obsm_pd'] = adata
