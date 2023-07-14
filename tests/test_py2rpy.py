from __future__ import annotations

from typing import TYPE_CHECKING
from warnings import catch_warnings, simplefilter

import numpy as np
import pytest
import scanpy as sc
from anndata import AnnData
from pandas import DataFrame
from rpy2.robjects import baseenv, globalenv
from rpy2.robjects.conversion import localconverter

import anndata2ri
from anndata2ri._rpy2_ext import importr


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp

    from anndata2ri.test_utils import Py2R


def mk_ad_simple() -> AnnData:
    return AnnData(
        np.array([[1, 2, 3], [0.3, 0.2, 0.1]]),
        dict(cluster=[1, 2]),
        obsm=dict(X_pca=np.array([[1, 3, 5, 7], [2, 4, 6, 8]])),
    )


def check_empty(_: Sexp) -> None:
    pass


def check_pca(ex: Sexp) -> None:
    sce = importr('SingleCellExperiment')
    assert [str(n) for n in sce.reducedDimNames(ex)] == ['PCA']
    pca = sce.reducedDim(ex, 'PCA')
    assert tuple(baseenv['dim'](pca)) == (2, 4)


datasets = [
    pytest.param(check_empty, (0, 0), AnnData, id='empty'),
    pytest.param(check_pca, (2, 3), mk_ad_simple, id='simple'),
    pytest.param(check_empty, (640, 11), sc.datasets.krumsiek11, id='krumsiek'),
]


@pytest.mark.parametrize(('check', 'shape', 'dataset'), datasets)
def test_py2rpy(
    py2r: Py2R,
    check: Callable[[Sexp], None],
    shape: tuple[int, ...],
    dataset: Callable[[], AnnData],
) -> None:
    if dataset is sc.datasets.krumsiek11:
        with pytest.warns(UserWarning, match=r'Duplicated obs_names'):
            ex = py2r(anndata2ri, dataset())
    else:
        ex = py2r(anndata2ri, dataset())
    assert tuple(baseenv['dim'](ex)[::-1]) == shape
    check(ex)


def test_py2rpy2_numpy_pbmc68k() -> None:
    """Not tested above as the pbmc68k dataset has some weird metadata."""
    from scanpy.datasets import pbmc68k_reduced

    try:
        anndata2ri.activate()
        with catch_warnings(record=True) as logs:
            simplefilter('ignore', DeprecationWarning)
            globalenv['adata'] = pbmc68k_reduced()
        assert len(logs) == 0, [m.message for m in logs]
    finally:
        anndata2ri.deactivate()


@pytest.mark.parametrize('attr', ['X', 'layers', 'obsm'])
def test_dfs(attr: str) -> None:
    """X, layers, obsm can contain dataframes."""
    adata = mk_ad_simple()
    if attr == 'X':
        adata.X = DataFrame(adata.X, index=adata.obs_names)
    elif attr == 'layers':
        adata.layers['X2'] = DataFrame(adata.X, index=adata.obs_names)
    elif attr == 'obsm':
        adata.obsm['X_pca'] = DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
    else:
        pytest.fail('Forgot to add a case')

    with localconverter(anndata2ri.converter):
        globalenv['adata_obsm_pd'] = adata


def test_df_error() -> None:
    adata = mk_ad_simple()
    adata.obsm['stuff'] = DataFrame(dict(a=[1, 2], b=list('ab'), c=[1.0, 2.0]), index=adata.obs_names)
    with pytest.raises(ValueError, match=r"DataFrame contains non-numeric columns \['b'\]"):
        anndata2ri.converter.py2rpy(adata)
