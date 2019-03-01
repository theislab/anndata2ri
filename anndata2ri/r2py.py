from typing import Optional, Union

import pandas as pd
from anndata import AnnData

from rpy2.rinterface import NULLType, SexpS4
from rpy2.robjects import numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.robject import RSlots
from rpy2.robjects.packages import importr

from .conv import converter, full_converter


@converter.rpy2py.register(SexpS4)
def rpy2py_s4(obj: SexpS4) -> Optional[Union[pd.DataFrame, AnnData]]:
    """
    See here for the slots: https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
    """
    if "DataFrame" in obj.rclass:
        return rpy2py_data_frame(obj)
    elif "SingleCellExperiment" in obj.rclass:
        return rpy2py_single_cell_experiment(obj)
    else:
        return default_converter.rpy2py(obj)


def rpy2py_data_frame(obj: SexpS4) -> pd.DataFrame:
    """
    S4 DataFrame class, not data.frame
    """
    with localconverter(default_converter):
        slots = RSlots(obj)
        columns = dict(slots["listData"].items())
        rownames = slots["rownames"]
        if isinstance(rownames, NULLType):
            rownames = pd.RangeIndex(slots["nrows"][0])

    return pd.DataFrame(columns, index=rownames)


def rpy2py_single_cell_experiment(obj: SexpS4) -> AnnData:
    with localconverter(default_converter):
        s4v = importr("S4Vectors")
        se = importr("SummarizedExperiment")
        sce = importr("SingleCellExperiment")

        assay_names = se.assayNames(obj)
        if not isinstance(assay_names, NULLType):
            assay_names = [str(a) for a in se.assayNames(obj)]
            # The assays can be stored in an env or elsewise so we don’t use obj.slots['assays']
            assays = [numpy2ri.rpy2py(assay).T for assay in (se.assay(obj, n) for n in assay_names)]
            # There’s SingleCellExperiment with no assays
            exprs, layers = assays[0], dict(zip(assay_names[1:], assays[1:]))
            assert len(exprs.shape) == 2, exprs.shape
        else:
            exprs, layers = None, {}

        rdim_names = [str(t) for t in sce.reducedDimNames(obj)]
        if rdim_names:
            reduced_dims = [numpy2ri.rpy2py(rd) for rd in (sce.reducedDim(obj, t) for t in rdim_names)]
            obsm = dict(zip(rdim_names, reduced_dims))
        else:
            obsm = None

        col_data = se.colData(obj)
        row_data = se.rowData(obj)
        metadata = s4v.metadata(obj)

    obs = rpy2py_data_frame(col_data)
    var = rpy2py_data_frame(row_data)
    # The whole shebang: configured converter, numpy, pandas and ours
    with localconverter(full_converter()):
        uns = dict(metadata.items())

    return AnnData(exprs, obs, var, uns, obsm, layers=layers)
