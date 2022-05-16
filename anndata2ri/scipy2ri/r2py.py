from warnings import warn

import numpy as np
from rpy2.rinterface import SexpS4, NULL
from rpy2.robjects import default_converter, numpy2ri, baseenv
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.robject import RSlots
from scipy import sparse

from .support import supported_r_matrix_classes
from .conv import converter


@converter.rpy2py.register(SexpS4)
def rmat_to_spmat(rmat: SexpS4):
    slots = RSlots(rmat)
    with localconverter(default_converter + numpy2ri.converter):
        shape = baseenv["dim"](rmat)
        shape = None if shape is NULL else tuple(shape)
        r_classes = set(rmat.rclass)
        if not supported_r_matrix_classes() & r_classes:
            if any(c.endswith("Matrix") for c in r_classes):
                warn(f"Encountered Matrix class that is not supported: {r_classes}")
            return rmat
        for storage, mat_cls, idx, nnz in [
            ("C", sparse.csc_matrix, lambda: [slots["i"], slots["p"]], lambda c: len(c[0])),
            ("R", sparse.csr_matrix, lambda: [slots["j"], slots["p"]], lambda c: len(c[0])),
            ("T", sparse.coo_matrix, lambda: [(slots["i"], slots["j"])], lambda c: len(c[0][0])),
            ("di", sparse.dia_matrix, lambda: [[0]], None),
        ]:
            if not supported_r_matrix_classes(storage=storage) & r_classes:
                continue
            coord_spec = idx()
            if supported_r_matrix_classes(types="n") & r_classes:
                # we have pattern matrix without data (but always i and j!)
                data = np.repeat(True, nnz(coord_spec))
            else:
                data = slots["x"]
            return mat_cls((data, *coord_spec), shape=shape)
        else:
            assert False, "Should have hit one of the branches"
