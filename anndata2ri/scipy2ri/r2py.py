from warnings import warn

from rpy2.rinterface import SexpS4
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
        shape = tuple(baseenv["dim"](rmat))
        r_classes = set(rmat.rclass)
        for storage, mat_cls, idx in [
            ("C", sparse.csc_matrix, lambda: [slots["i"], slots["p"]]),
            ("R", sparse.csr_matrix, lambda: [slots["j"], slots["p"]]),
            ("T", sparse.coo_matrix, lambda: [(slots["i"], slots["j"])]),
            ("di", sparse.dia_matrix, lambda: [[0]]),
        ]:
            if supported_r_matrix_classes(storage=storage) & r_classes:
                return mat_cls((slots["x"], *idx()), shape=shape)
    if any(c.endswith("Matrix") for c in r_classes):
        warn(f"Encountered Matrix class that is not supported: {r_classes}")
    return rmat
