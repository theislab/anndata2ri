from collections.abc import Mapping
from warnings import warn

import numpy as np
import pandas as pd
from anndata import AnnData

from rpy2.robjects import conversion, pandas2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import ListVector
from rpy2.robjects.methods import RS4

from . import conv_name
from .conv import converter, mat_converter, full_converter
from .rpy2_ext import importr


class NotConvertedWarning(Warning):
    pass


dict_converter = conversion.Converter("Converter handling dicts")
dict_converter.py2rpy.register(np.bool_, lambda x: conversion.py2rpy(bool(x)))
dict_converter.py2rpy.register(np.int_, lambda x: conversion.py2rpy(int(x)))
dict_converter.py2rpy.register(np.float_, lambda x: conversion.py2rpy(float(x)))
dict_converter.py2rpy.register(np.bytes_, lambda x: conversion.py2rpy(bytes(x)))
dict_converter.py2rpy.register(np.str_, lambda x: conversion.py2rpy(str(x)))


@dict_converter.py2rpy.register(Mapping)
def py2rpy_dict(obj: Mapping) -> ListVector:
    """Try converting everything. For nested dicts, this needs itself to be registered"""
    converted = {}
    for k, v in obj.items():
        try:
            converted[str(k)] = conversion.py2rpy(v)
        except NotImplementedError as e:
            warn(str(e), NotConvertedWarning)
    # This tries to convert everything again. This works because py2rpy(Sexp) is the identity function
    return ListVector(converted)


def check_no_dupes(idx: pd.Index, name: str) -> bool:
    dupes = idx.duplicated().any()
    if dupes:
        warn(f"Duplicated {name}: {idx[idx.duplicated(False)].sort_values()}")
    return not dupes


@converter.py2rpy.register(AnnData)
def py2rpy_anndata(obj: AnnData) -> RS4:
    with localconverter(default_converter):
        s4v = importr("S4Vectors")
        sce = importr("SingleCellExperiment")
        # TODO: sparse
        x = {} if obj.X is None else dict(X=mat_converter.py2rpy(obj.X.T))
        layers = {k: mat_converter.py2rpy(v.T) for k, v in obj.layers.items()}
        assays = ListVector({**x, **layers})

        row_args = {k: pandas2ri.py2rpy(v) for k, v in obj.var.items()}
        if check_no_dupes(obj.var_names, "var_names"):
            row_args["row.names"] = pandas2ri.py2rpy(obj.var_names)
        row_data = s4v.DataFrame(**row_args)

        col_args = {k: pandas2ri.py2rpy(v) for k, v in obj.obs.items()}
        if check_no_dupes(obj.obs_names, "obs_names"):
            col_args["row.names"] = pandas2ri.py2rpy(obj.obs_names)
        col_data = s4v.DataFrame(**col_args)

        # Convert everything we know
        with localconverter(full_converter() + dict_converter):
            metadata = ListVector(obj.uns.items())

        rd_args = {conv_name.scanpy2sce(k): mat_converter.py2rpy(obj.obsm[k]) for k in obj.obsm.keys()}
        reduced_dims = s4v.SimpleList(**rd_args)

        return sce.SingleCellExperiment(
            assays=assays, rowData=row_data, colData=col_data, metadata=metadata, reducedDims=reduced_dims
        )
