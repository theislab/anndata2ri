from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from rpy2.rinterface import NULL, Sexp, RTYPES
from rpy2.robjects import conversion, default_converter, numpy2ri, pandas2ri
from rpy2.robjects.vectors import Matrix, ListVector
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr


converter = conversion.Converter("original anndata conversion")
rpy2py = converter.rpy2py
py2rpy = converter.py2rpy


# Python to R


@py2rpy.register(AnnData)
def py2rpy_anndata(obj: AnnData) -> RS4:
    s4v = importr("S4Vectors")
    sce = importr("SingleCellExperiment")

    layers = {k: conversion.py2rpy(v.T) for k, v in obj.layers.items()}
    assays = ListVector({"X": conversion.py2rpy(obj.X.T), **layers})

    row_args = {k: conversion.py2rpy(v) for k, v in obj.var.items()}
    row_args["row.names"] = conversion.py2rpy(obj.var.index)
    row_data = s4v.DataFrame(**row_args)

    col_args = {k: conversion.py2rpy(v) for k, v in obj.obs.items()}
    col_args["row.names"] = conversion.py2rpy(obj.obs.index)
    col_data = s4v.DataFrame(**col_args)

    metadata = ListVector({k: conversion.py2rpy(v) for k, v in obj.uns.items()})

    return sce.SingleCellExperiment(assays=assays, rowData=row_data, colData=col_data, metadata=metadata)


# R to Python


# https://bitbucket.org/rpy2/rpy2/issues/518/converting-matrices-is-very-slow
@numpy2ri.rpy2py.register(Sexp)  # No idea why this is necessary. I’d think the dispatcher catches this earlier
@rpy2py.register(Matrix)
def rpy2py_matrix_ad(obj: Sexp):
    """
    For some reason the original matrix conversion is dog slow.
    Using memoryview fixes that.
    """
    if obj.typeof in numpy2ri._vectortypes and obj.typeof != RTYPES.VECSXP:
        if hasattr(obj, "memoryview"):
            return np.asarray(obj.memoryview())
        else:
            return np.asarray(obj)
    else:  # If i used the pandas converted, this’d be a recursion error
        return default_converter.rpy2py(obj)


@rpy2py.register(RS4)
def rpy2py_s4(obj: RS4) -> Optional[Union[pd.DataFrame, AnnData]]:
    """
    See here for the slots: https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
    """
    if "DataFrame" in obj.rclass:
        return rpy2py_data_frame(obj)
    elif "SingleCellExperiment" in obj.rclass:
        return rpy2py_single_cell_experiment(obj)
    # else default to None


def rpy2py_data_frame(obj: RS4) -> pd.DataFrame:
    """
    S4 DataFrame class, not data.frame
    """
    columns = {k: rpy2py(v) if isinstance(v, Sexp) else v for k, v in obj.slots["listData"].items()}
    rownames = obj.slots["rownames"]
    if rownames is NULL:
        rownames = None

    return pd.DataFrame(columns, index=rownames)


def rpy2py_single_cell_experiment(obj: RS4) -> AnnData:
    se = importr("SummarizedExperiment")
    # sce = importr('SingleCellExperiment')

    assay_names = se.assayNames(obj)
    if assay_names is not NULL:
        # The assays can be stored in an env or elsewise so we don’t use obj.slots['assays']
        assays = [se.assay(obj, str(n)).T for n in assay_names]
        # There’s SingleCellExperiment with no assays
        exprs, layers = assays[0], dict(zip(assay_names[1:], assays[1:]))
        assert len(exprs.shape) == 2, exprs.shape
    else:
        exprs, layers = None, {}

    obs = rpy2py(se.colData(obj))
    assert isinstance(obs, pd.DataFrame), type(obs)
    var = rpy2py(se.rowData(obj))
    assert isinstance(var, pd.DataFrame), type(var)

    # TODO: se.metadata, se.dimnames

    return AnnData(exprs, obs, var, layers=layers)


# Activation / deactivation


def create_converter() -> conversion.Converter:
    pandas2ri.activate()
    new_converter = conversion.Converter("anndata conversion", template=conversion.converter)
    pandas2ri.deactivate()

    for k, v in py2rpy.registry.items():
        if k is not object:
            new_converter.py2rpy.register(k, v)

    for k, v in rpy2py.registry.items():
        if k is not object:
            new_converter.rpy2py.register(k, v)

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
