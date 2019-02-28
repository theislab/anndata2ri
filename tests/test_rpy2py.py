import pytest
import pandas as pd
from anndata import AnnData
from rpy2.robjects import r, default_converter
from rpy2.robjects.conversion import ConversionContext
from rpy2.robjects.packages import importr, data

import anndata2ri


se = importr("SummarizedExperiment")
sce = importr("SingleCellExperiment")
sc_rna_seq_data = data(importr("scRNAseq"))
as_ = getattr(importr("methods"), "as")


sumex_allen = sc_rna_seq_data.fetch("allen")["allen"]
ex_allen = (379, 20908), lambda: as_(sumex_allen, "SingleCellExperiment")
ex_empty = (0, 0), lambda: sce.SingleCellExperiment()


@pytest.mark.parametrize("shape,dataset", [ex_empty, ex_allen])
def test_convert_manual(shape, dataset):
    ad = anndata2ri.converter.rpy2py(dataset())
    assert isinstance(ad, AnnData)
    assert ad.shape == shape


@pytest.mark.parametrize("shape,dataset", [ex_empty, ex_allen])
def test_convert_with(shape, dataset):
    # Needs default_converter to call `as` on the SummarizedExperiment:
    # Calling a R function returning a S4 object requires py2rpy[RS4], py2rpy[str], …
    with ConversionContext(default_converter + anndata2ri.converter):
        ad = dataset()
    assert isinstance(ad, AnnData)
    assert ad.shape == shape


@pytest.mark.parametrize("shape,dataset", [ex_empty, ex_allen])
def test_convert_activate(shape, dataset):
    try:
        anndata2ri.activate()
        ad = dataset()
    finally:
        anndata2ri.deactivate()
    assert isinstance(ad, AnnData)
    assert ad.shape == shape


def test_convert_empty_df_with_rows():
    df = r("S4Vectors::DataFrame(a=1:10)[, -1]")
    assert df.slots["nrows"][0] == 10

    c = anndata2ri.create_converter()
    df_py = c.rpy2py(df)
    assert isinstance(df_py, pd.DataFrame)
