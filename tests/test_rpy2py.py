from pathlib import Path

import pandas as pd
import pytest
from anndata import AnnData
from rpy2.robjects import conversion, r

import anndata2ri
from anndata2ri.rpy2_ext import importr


as_ = getattr(importr('methods'), 'as')
se = importr('SummarizedExperiment')
sce = importr('SingleCellExperiment')
eh = importr('ExperimentHub')
seq = importr('scRNAseq')


# avoid prompt
Path(eh.getExperimentHubOption('CACHE')[0]).mkdir(parents=True, exist_ok=True)


def check_allen(adata):
    assert adata.uns.keys() == {'SuppInfo', 'which_qc'}
    assert set(adata.obs.keys()) > {'NREADS', 'NALIGNED', 'Animal.ID', 'passes_qc_checks_s'}
    assert adata.obs['Secondary.Type'][:4].tolist() == ['L4 Ctxn3', '', 'L5a Batf3', None], 'NAs not conserved?'
    assert adata.obs['Animal.ID'][:4].tolist() == [133632, 133632, 151560, pd.NA], 'NAs not conserved?'


def check_example(adata):
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
    pytest.param(
        check_allen,
        (379, 20816),
        lambda: as_(seq.ReprocessedAllenData(assays='tophat_counts'), 'SingleCellExperiment'),
        id='allen',
    ),
    pytest.param(lambda x: None, (0, 0), sce.SingleCellExperiment, id='empty'),
    pytest.param(check_example, (100, 200), lambda: r(code_example), id='example'),
]


@pytest.mark.parametrize('check,shape,dataset', expression_sets)
def test_convert_manual(r2py, check, shape, dataset):
    ad = r2py(anndata2ri, dataset)
    assert isinstance(ad, AnnData)
    assert ad.shape == shape
    check(ad)


def test_convert_empty_df_with_rows(r2py):
    df = r('S4Vectors::DataFrame(a=1:10)[, -1]')
    assert df.slots['nrows'][0] == 10

    df_py = r2py(anndata2ri, lambda: conversion.rpy2py(df))
    assert isinstance(df_py, pd.DataFrame)


def test_convert_factor(r2py):
    code = """
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(rpois(6*4, 5), ncol=4)),
        colData = S4Vectors::DataFrame(a_factor = factor(c(rep('A', 2), NA, rep('B', 1))))
    )
    """
    ad = r2py(anndata2ri, lambda: r(code))
    assert isinstance(ad.obs['a_factor'].values, pd.Categorical)
    assert ad.obs['a_factor'].values.tolist() == pd.Categorical.from_codes([0, 0, -1, 1], ['A', 'B']).tolist()
