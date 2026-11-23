from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import h5py
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from anndata.abc import CSCDataset, CSRDataset
from rpy2.robjects import baseenv, conversion, numpy2ri, r
from scipy import sparse

import anndata2ri
from anndata2ri._rpy2_ext import importr


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp

    from anndata2ri.test_utils import R2Py


as_ = getattr(importr('methods'), 'as')
se = importr('SummarizedExperiment')
sce = importr('SingleCellExperiment')


def get_allen() -> Sexp:
    eh = importr('ExperimentHub')
    seq = importr('scRNAseq')

    # avoid prompt
    Path(eh.getExperimentHubOption('CACHE')[0]).mkdir(parents=True, exist_ok=True)
    return as_(seq.ReprocessedAllenData(assays='tophat_counts'), 'SingleCellExperiment')


def check_allen(adata: AnnData) -> None:
    assert adata.uns.keys() == {'SuppInfo', 'which_qc'}
    assert set(adata.obs.keys()) > {'NREADS', 'NALIGNED', 'Animal.ID', 'passes_qc_checks_s'}
    sec_type = adata.obs['Secondary.Type'][:4].tolist()
    assert sec_type[:3] == ['L4 Ctxn3', '', 'L5a Batf3'], 'NAs not conserved?'
    assert pd.isna(sec_type[3]), 'NAs not conserved?'
    assert adata.obs['Animal.ID'][:4].tolist() == [133632, 133632, 151560, pd.NA], 'NAs not conserved?'


def check_example(adata: AnnData) -> None:
    assert set(adata.obsm.keys()) == {'X_pca', 'X_tsne'}
    assert adata.obsm['X_pca'].shape == (100, 5)


code_example = """
local({
    ncells <- 100
    u <- matrix(rpois(20000, 5), ncol=ncells)
    p <- matrix(runif(ncells*5), ncells)
    t <- matrix(rnorm(ncells*2), ncells)
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = u, logcounts = log2(u + 1)),
        reducedDims = S4Vectors::SimpleList(PCA = p, tSNE = t)
    )
})
"""

expression_sets = [
    pytest.param(check_allen, (379, 20816), get_allen, id='allen'),
    pytest.param(lambda _: None, (0, 0), sce.SingleCellExperiment, id='empty'),
    pytest.param(check_example, (100, 200), lambda: r(code_example), id='example'),
]


@pytest.mark.parametrize(('check', 'shape', 'dataset'), expression_sets)
def test_manual(
    r2py: R2Py,
    check: Callable[[AnnData], None],
    shape: tuple[int, ...],
    dataset: Callable[[], Sexp],
) -> None:
    ad = r2py(anndata2ri, dataset)
    assert isinstance(ad, AnnData)
    assert ad.shape == shape
    check(ad)


def test_empty_df_with_rows(r2py: R2Py) -> None:
    df = r('S4Vectors::DataFrame(a=1:10)[, -1]')
    assert df.slots['nrows'][0] == 10

    df_py = r2py(anndata2ri, lambda: conversion.get_conversion().rpy2py(df))
    assert isinstance(df_py, pd.DataFrame)


def test_factor(r2py: R2Py) -> None:
    code = """
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(rpois(6*4, 5), ncol=4)),
        colData = S4Vectors::DataFrame(a_factor = factor(c(rep('A', 2), NA, rep('B', 1))))
    )
    """
    ad = r2py(anndata2ri, lambda: r(code))
    assert isinstance(ad.obs['a_factor'].values, pd.Categorical)
    assert ad.obs['a_factor'].to_numpy().tolist() == pd.Categorical.from_codes([0, 0, -1, 1], ['A', 'B']).tolist()


delayed_mats = [
    pytest.param(np.ndarray, 'DelayedArray::DelayedArray(matrix((1:6) - 3, 2))', id='dense'),
    pytest.param(
        sparse.csc_array,
        'DelayedArray::DelayedArray(Matrix::sparseMatrix(1:2, 3:2, x = c(1, -2), dims = c(4, 3)))',
        id='sparse',
    ),
    pytest.param(
        np.ndarray,
        'BiocSingular::LowRankMatrix(matrix((1:6) / 2, 3), matrix((1:8) - 4, 4))',
        id='LowRankMatrix',
    ),
]


@pytest.mark.parametrize(('cls', 'code'), delayed_mats)
def test_delayed(r2py: R2Py, cls: type, code: str) -> None:
    converted = r2py(anndata2ri, lambda: r(code))
    assert isinstance(converted, cls)
    dense = numpy2ri.converter.rpy2py(baseenv['as.matrix'](r(code)))
    np.testing.assert_allclose(converted.toarray() if isinstance(converted, sparse.sparray) else converted, dense)


def test_delayed_assays(r2py: R2Py) -> None:
    code = """
    local({
        lr <- BiocSingular::LowRankMatrix(matrix(rnorm(12), 6), matrix(rnorm(8), 4))
        SingleCellExperiment::SingleCellExperiment(
            assays = list(reconstructed = lr, delayed = DelayedArray::DelayedArray(as.matrix(lr)))
        )
    })
    """
    ad = r2py(anndata2ri, lambda: r(code))
    assert ad.shape == (4, 6)
    assert isinstance(ad.X, np.ndarray)
    np.testing.assert_allclose(ad.layers['delayed'], ad.X)


X_BACKED = np.array([[0.0, 1.0], [0.0, 2.0], [3.0, 0.0]])  # cells × genes
R_X_BACKED = 'matrix(c(0, 0, 3, 1, 2, 0), 3)'


def write_r(code: str) -> Callable[[Path], object]:
    return lambda p: r(code.format(x=R_X_BACKED, path=p.as_posix()))


def write_h5ad(fmt: type[sparse.sparray] | None) -> Callable[[Path], None]:
    return lambda p: AnnData(X_BACKED if fmt is None else fmt(X_BACKED)).write_h5ad(p)


backed_assays = [
    pytest.param(
        write_r('HDF5Array::writeHDF5Array(t({x}), "{path}", "m")'),
        'HDF5Array::HDF5Array("{path}", "m")',
        h5py.Dataset,
        id='HDF5Array',
    ),
    pytest.param(
        write_r('HDF5Array::writeTENxMatrix(Matrix::Matrix(t({x}), sparse = TRUE), "{path}", "g")'),
        'HDF5Array::TENxMatrix("{path}", "g")',
        CSRDataset,
        id='TENxMatrix',
    ),
    pytest.param(write_h5ad(None), 'HDF5Array::H5ADMatrix("{path}")', h5py.Dataset, id='h5ad-dense'),
    pytest.param(write_h5ad(sparse.csr_array), 'HDF5Array::H5ADMatrix("{path}")', CSRDataset, id='h5ad-csr'),
    pytest.param(write_h5ad(sparse.csc_array), 'HDF5Array::H5ADMatrix("{path}")', CSCDataset, id='h5ad-csc'),
    pytest.param(
        write_r('HDF5Array::writeHDF5Array(t({x}), "{path}", "m")'),
        '1 * HDF5Array::HDF5Array("{path}", "m")',
        np.ndarray,
        id='not-pristine',
    ),
]


@pytest.mark.parametrize(('write', 'read', 'cls'), backed_assays)
def test_backed_assays(r2py: R2Py, tmp_path: Path, write: Callable[[Path], object], read: str, cls: type) -> None:
    path = tmp_path / 'x.h5'
    write(path)
    read = read.format(path=path.as_posix())
    # reading the same file twice must work too
    code = f'SingleCellExperiment::SingleCellExperiment(list(counts = {read}, again = {read}))'
    ad = r2py(anndata2ri, lambda: r(code))
    for arr in [ad.X, ad.layers['again']]:
        assert isinstance(arr, cls)
        x = arr[()] if isinstance(arr, h5py.Dataset) else arr
        x = x.to_memory() if isinstance(x, CSRDataset | CSCDataset) else x
        np.testing.assert_array_equal(x.toarray() if sparse.issparse(x) else x, X_BACKED)
