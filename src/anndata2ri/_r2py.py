from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pandas as pd
from anndata import AnnData
from rpy2.rinterface import IntSexpVector, NULLType, Sexp, SexpS4, baseenv
from rpy2.robjects import RS4, default_converter, numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.robject import RSlots

from . import _conv_name
from ._conv import converter, mat_rpy2py
from ._rpy2_ext import R_INT_BYTES, importr
from .scipy2ri import supported_r_matrix_classes
from .scipy2ri._r2py import rmat_to_spmat


if TYPE_CHECKING:
    from collections.abc import Mapping

    import numpy as np
    from scipy.sparse import spmatrix


@converter.rpy2py.register(SexpS4)
def rpy2py_s4(obj: SexpS4) -> pd.DataFrame | AnnData | None:
    """Convert known S4 class instance to Python object.

    See here for the slots: https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
    """
    r_classes = set(obj.rclass)
    if {'DataFrame', 'DFrame'} & r_classes:
        return rpy2py_data_frame(obj)
    if 'SingleCellExperiment' in r_classes:
        return rpy2py_single_cell_experiment(obj)
    if supported_r_matrix_classes() & r_classes:
        return rmat_to_spmat(obj)
    # Don’t use the registered one, it would lead to recursion.
    return default_converter.rpy2py(obj)


def rpy2py_vector(v: Sexp) -> Any:  # noqa: ANN401
    """Convert R vector to Python vectors.

    Also handles NA in int vectors: https://github.com/rpy2/rpy2/issues/376
    """
    if not isinstance(v, Sexp):
        return v
    if isinstance(v, IntSexpVector):
        assert v._R_SIZEOF_ELT == R_INT_BYTES, 'R integer size changed away from 32 bit'  # noqa: SLF001
        r = pd.array(v, dtype=pd.Int32Dtype())
        v_is_na = numpy2ri.rpy2py(baseenv['is.na'](v)).astype(bool)
        if 'factor' in v.rclass:
            levels = numpy2ri.rpy2py(baseenv['levels'](v))
            codes = r.to_numpy() - 1
            # temporarily set NA values to a valid index
            codes[v_is_na] = 0
            codes = codes.astype(int)
            r = pd.array(levels[codes], dtype=pd.CategoricalDtype(levels))
        r[v_is_na] = pd.NA
        return r
    return pandas2ri.rpy2py(v)


def rpy2py_data_frame(obj: SexpS4) -> pd.DataFrame:
    """S4 DataFrame class, not data.frame."""
    slots = RSlots(obj)
    with localconverter(default_converter):
        columns = {k: rpy2py_vector(v) for k, v in slots['listData'].items()}
        rownames = slots['rownames']
        if isinstance(rownames, NULLType):
            rownames = pd.RangeIndex(slots['nrows'][0])

    return pd.DataFrame(columns, index=rownames)


def rpy2py_single_cell_experiment(obj: SexpS4) -> AnnData:
    with localconverter(default_converter):
        s4v = importr('S4Vectors')
        se = importr('SummarizedExperiment')
        sce = importr('SingleCellExperiment')

        def convert_mats(
            attr: str, mats: Mapping[str, Sexp], *, transpose: bool = False
        ) -> list[np.ndarray | spmatrix]:
            rv = []
            for n, mat in mats.items():
                conv = mat_rpy2py(mat)
                if isinstance(conv, RS4):
                    cls_names = mat_rpy2py(conv.slots['class']).tolist()
                    msg = f'Cannot convert {attr} “{n}” of type(s) {cls_names} to Python'
                    raise TypeError(msg)
                rv.append(conv.T if transpose else conv)
            return rv

        assay_names = se.assayNames(obj)
        if not isinstance(assay_names, NULLType):
            assay_names = [str(a) for a in se.assayNames(obj)]
            # The assays can be stored in an env or elsewise so we don’t use obj.slots['assays']
            assays = convert_mats('assay', {n: se.assay(obj, n) for n in assay_names}, transpose=True)
            # There’s SingleCellExperiment with no assays
            exprs, layers = assays[0], dict(zip(assay_names[1:], assays[1:]))
            assert len(exprs.shape) == 2, exprs.shape  # noqa: PLR2004
        else:
            exprs, layers = None, {}

        rdim_names = sce.reducedDimNames(obj)
        if not isinstance(rdim_names, NULLType):
            rdim_names = [str(t) for t in rdim_names]
            reduced_dims = convert_mats('reducedDim', {t: sce.reducedDim(obj, t) for t in rdim_names})
            obsm = {_conv_name.sce2scanpy(n): d for n, d in zip(rdim_names, reduced_dims)}
        else:
            obsm = None

        col_data = se.colData(obj)
        row_data = se.rowData(obj)
        metadata = s4v.metadata(obj)

    obs = rpy2py_data_frame(col_data)
    var = rpy2py_data_frame(row_data)
    # To avoid ImplicitModificationWarning
    obs.index = obs.index.astype('string')
    var.index = var.index.astype('string')
    # The whole shebang: configured converter, numpy, pandas and ours
    with localconverter(converter):
        uns = dict(metadata.items())

    return AnnData(exprs, obs, var, uns, obsm, layers=layers)
