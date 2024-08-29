from __future__ import annotations

from typing import TYPE_CHECKING
from warnings import catch_warnings, filterwarnings, simplefilter

import numpy as np
import pytest
import scanpy as sc
from anndata import AnnData
from pandas import DataFrame
from rpy2.robjects import baseenv, globalenv
from rpy2.robjects.conversion import localconverter
from scipy import sparse

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


@pytest.mark.parametrize('dtype', [np.float32, np.float64, np.int32, np.int64])
@pytest.mark.parametrize('mat_type', [np.asarray, sparse.csr_matrix])
def test_simple(
    py2r: Py2R,
    dtype: np.dtype,
    mat_type: Callable[[np.ndarray, np.dtype], np.ndarray | sparse.spmatrix],
) -> None:
    data = mk_ad_simple()
    if data.X is not None:
        data.X = mat_type(data.X, dtype=dtype)
    ex = py2r(anndata2ri, data)
    assert tuple(baseenv['dim'](ex)[::-1]) == data.shape


def krumsiek() -> AnnData:
    adata = sc.datasets.krumsiek11()
    adata.obs_names_make_unique()
    return adata


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
    pytest.param(check_empty, (640, 11), krumsiek, id='krumsiek'),
]


@pytest.mark.parametrize(('check', 'shape', 'dataset'), datasets)
def test_datasets(
    py2r: Py2R,
    check: Callable[[Sexp], None],
    shape: tuple[int, ...],
    dataset: Callable[[], AnnData],
) -> None:
    if dataset is krumsiek:
        # TODO(flying-sheep): Adapt to rpy2 changes instead
        # https://github.com/theislab/anndata2ri/issues/109
        with pytest.warns(DeprecationWarning, match=r'rpy2\.robjects\.conversion is deprecated'):
            filterwarnings('ignore', r'Duplicated obs_names', UserWarning)
            filterwarnings('ignore', r'Observation names are not unique', UserWarning)
            ex = py2r(anndata2ri, dataset())
    else:
        ex = py2r(anndata2ri, dataset())
    assert tuple(baseenv['dim'](ex)[::-1]) == shape
    check(ex)


def test_numpy_pbmc68k() -> None:
    """Not tested above as the pbmc68k dataset has some weird metadata."""
    from scanpy.datasets import pbmc68k_reduced

    try:
        with pytest.warns(DeprecationWarning, match=r'global conversion'):
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
