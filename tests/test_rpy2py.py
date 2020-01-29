import pytest
import pandas as pd
from anndata import AnnData
from rpy2.robjects import r, conversion

import anndata2ri
from anndata2ri.rpy2_ext import importr, data
from anndata2ri.test_utils import conversions_rpy2py

as_ = getattr(importr("methods"), "as")
se = importr("SummarizedExperiment")
sce = importr("SingleCellExperiment")
sumex_allen = data("scRNAseq", "allen")["allen"]


def check_allen(adata):
    assert adata.uns.keys() == {"SuppInfo", "which_qc"}
    assert set(adata.obs.keys()) > {"NREADS", "NALIGNED", "Animal.ID", "passes_qc_checks_s"}


def check_example(adata):
    assert set(adata.obsm.keys()) == {"X_pca", "X_tsne"}
    assert adata.obsm["X_pca"].shape == (100, 5)


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
    pytest.param(check_allen, (379, 20908), lambda: as_(sumex_allen, "SingleCellExperiment"), id="allen"),
    pytest.param(lambda x: None, (0, 0), sce.SingleCellExperiment, id="empty"),
    pytest.param(check_example, (100, 200), lambda: r(code_example), id="example"),
]


@pytest.mark.parametrize("convert", conversions_rpy2py)
@pytest.mark.parametrize("check,shape,dataset", expression_sets)
def test_convert_manual(convert, check, shape, dataset):
    ad = convert(anndata2ri, dataset)
    assert isinstance(ad, AnnData)
    assert ad.shape == shape
    check(ad)


@pytest.mark.parametrize("convert", conversions_rpy2py)
def test_convert_empty_df_with_rows(convert):
    df = r("S4Vectors::DataFrame(a=1:10)[, -1]")
    assert df.slots["nrows"][0] == 10

    df_py = convert(anndata2ri, lambda: conversion.rpy2py(df))
    assert isinstance(df_py, pd.DataFrame)


@pytest.mark.parametrize("convert", conversions_rpy2py)
def test_convert_factor(convert):
    code = """
    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(rpois(6*4, 5), ncol=4)),
        colData = S4Vectors::DataFrame(a_factor = factor(c(rep('A', 3), rep('B', 1))))
    )
    """
    ad = convert(anndata2ri, lambda: r(code))
    assert isinstance(ad.obs["a_factor"].values, pd.Categorical)
    assert all(ad.obs["a_factor"].values == pd.Categorical.from_codes([0, 0, 0, 1], ["A", "B"]))
