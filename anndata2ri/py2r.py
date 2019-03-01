from warnings import warn

import pandas as pd
from anndata import AnnData

from rpy2.robjects import conversion, pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import ListVector
from rpy2.robjects.methods import RS4
from rpy2.robjects.packages import importr

from .conv import converter, full_converter


dict_converter = conversion.Converter("Converter handling dicts")


@dict_converter.py2rpy.register(dict)
def py2rpy_dict(obj):
    return ListVector({str(k): v for k, v in obj.items()})


def check_no_dupes(idx: pd.Index, name: str):
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
        x = {} if obj.X is None else dict(X=numpy2ri.py2rpy(obj.X.T))
        layers = {k: numpy2ri.py2rpy(v.T) for k, v in obj.layers.items()}
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

        return sce.SingleCellExperiment(assays=assays, rowData=row_data, colData=col_data, metadata=metadata)
