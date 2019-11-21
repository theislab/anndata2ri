import pytest
import pandas as pd
from anndata import AnnData
from rpy2.robjects import r, conversion
from rpy2.robjects.packages import importr, data

import anndata2ri
from anndata2ri.test_utils import conversions_rpy2py

se = importr("SummarizedExperiment")
sce = importr("SingleCellExperiment")
sc_rna_seq_data = data(importr("scRNAseq"))
as_ = getattr(importr("methods"), "as")


def check_allen(adata):
    assert adata.uns.keys() == {"SuppInfo", "which_qc"}
    assert set(adata.obs.keys()) > {"NREADS", "NALIGNED", "Animal.ID", "passes_qc_checks_s"}


def check_example(adata):
    assert set(adata.obsm.keys()) == {"X_pca", "X_tsne"}
    assert adata.obsm["X_pca"].shape == (100, 5)


sumex_allen = sc_rna_seq_data.fetch("allen")["allen"]
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
