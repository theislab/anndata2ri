import pytest
from rpy2.robjects import conversion
from rpy2.robjects.packages import importr, data

import anndata2ri


se = importr("SummarizedExperiment")
sce = importr("SingleCellExperiment")
sc_rna_seq_data = data(importr("scRNAseq"))
as_ = getattr(importr("methods"), "as")


ex_empty = sce.SingleCellExperiment()
ex_allen = as_(sc_rna_seq_data.fetch("allen")["allen"], "SingleCellExperiment")


@pytest.mark.parametrize("dataset", [ex_empty, ex_allen])
def test_convert_with(dataset):
    with conversion.ConversionContext(anndata2ri.create_converter()) as c:
        ad = c.rpy2py(dataset)
        print(ad)


@pytest.mark.parametrize("dataset", [ex_empty, ex_allen])
def test_convert_activate(dataset):
    try:
        anndata2ri.activate()
        ad = conversion.rpy2py(dataset)
        print(ad)
    finally:
        anndata2ri.deactivate()
