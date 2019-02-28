"""
Converter between Python’s AnnData and R’s SingleCellExperiment.
"""

from typing import Optional, Union

import pandas as pd
from anndata import AnnData
from get_version import get_version

from rpy2.rinterface import NULLType, SexpS4
from rpy2.robjects import conversion, pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter, overlay_converter
from rpy2.robjects.robject import RSlots
from rpy2.robjects.vectors import ListVector
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr


__version__ = get_version(__file__)

converter = conversion.Converter("original anndata conversion")
rpy2py = converter.rpy2py
py2rpy = converter.py2rpy


# Python to R


@py2rpy.register(AnnData)
def py2rpy_anndata(obj: AnnData) -> RS4:
    with localconverter(default_converter):
        s4v = importr("S4Vectors")
        sce = importr("SingleCellExperiment")
        # TODO: sparse
        layers = {k: numpy2ri.py2rpy(v.T) for k, v in obj.layers.items()}
        assays = ListVector({"X": numpy2ri.py2rpy(obj.X.T), **layers})

        row_args = {k: pandas2ri.py2rpy(v) for k, v in obj.var.items()}
        row_args["row.names"] = pandas2ri.py2rpy(obj.var_names)
        row_data = s4v.DataFrame(**row_args)

        col_args = {k: pandas2ri.py2rpy(v) for k, v in obj.obs.items()}
        col_args["row.names"] = pandas2ri.py2rpy(obj.obs_names)
        col_data = s4v.DataFrame(**col_args)

        metadata = ListVector({k: numpy2ri.py2rpy(v) for k, v in obj.uns.items()})

        return sce.SingleCellExperiment(assays=assays, rowData=row_data, colData=col_data, metadata=metadata)


# R to Python


@rpy2py.register(SexpS4)
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
        se = importr("SummarizedExperiment")
        # sce = importr('SingleCellExperiment')

        assay_names = se.assayNames(obj)
        if not isinstance(assay_names, NULLType):
            # The assays can be stored in an env or elsewise so we don’t use obj.slots['assays']
            assays = [numpy2ri.rpy2py(assay).T for assay in (se.assay(obj, str(n)) for n in assay_names)]
            # There’s SingleCellExperiment with no assays
            exprs, layers = assays[0], dict(zip(assay_names[1:], assays[1:]))
            assert len(exprs.shape) == 2, exprs.shape
        else:
            exprs, layers = None, {}

        col_data = se.colData(obj)
        row_data = se.rowData(obj)

    obs = rpy2py_data_frame(col_data)
    var = rpy2py_data_frame(row_data)

    # TODO: se.metadata

    return AnnData(exprs, obs, var, layers=layers)


# Activation / deactivation


def create_converter() -> conversion.Converter:
    pandas2ri.activate()
    new_converter = conversion.Converter("anndata conversion", template=conversion.converter)
    pandas2ri.deactivate()

    overlay_converter(converter, new_converter)

    return new_converter


original_converter: Optional[conversion.Converter] = None


def activate():
    global original_converter

    if original_converter is not None:
        return

    new_converter = create_converter()
    original_converter = conversion.converter
    conversion.set_conversion(new_converter)


def deactivate():
    global original_converter

    if original_converter is None:
        return

    conversion.set_conversion(original_converter)
    original_converter = None
