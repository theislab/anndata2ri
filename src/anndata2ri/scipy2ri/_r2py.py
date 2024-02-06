from __future__ import annotations

from typing import TYPE_CHECKING
from warnings import warn

import numpy as np
from rpy2.rinterface import NULL, SexpS4
from rpy2.robjects import baseenv, default_converter, numpy2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.robject import RSlots
from scipy import sparse

from ._conv import converter
from ._support import SupportedMatStor, supported_r_matrix_classes


if TYPE_CHECKING:
    from collections.abc import Callable, Iterable

    CoordSpec = tuple[np.ndarray, ...] | tuple[tuple[np.ndarray, ...]] | tuple[list[int]]


OPTIONS: Iterable[
    tuple[
        SupportedMatStor,
        type[sparse.spmatrix],
        Callable[[RSlots], CoordSpec],
        Callable[[CoordSpec], int] | None,
    ]
] = [
    ('C', sparse.csc_matrix, lambda slots: (slots['i'], slots['p']), lambda c: len(c[0])),
    ('R', sparse.csr_matrix, lambda slots: (slots['j'], slots['p']), lambda c: len(c[0])),
    ('T', sparse.coo_matrix, lambda slots: ((slots['i'], slots['j']),), lambda c: len(c[0][0])),
    ('di', sparse.dia_matrix, lambda _: ([0],), None),
]


@converter.rpy2py.register(SexpS4)
def rmat_to_spmat(rmat: SexpS4) -> sparse.spmatrix:
    """Convert R sparse matrices to scipy sparse matrices."""
    slots = RSlots(rmat)
    with localconverter(default_converter + numpy2ri.converter):
        shape = baseenv['dim'](rmat)
        shape = None if shape is NULL else tuple(shape)
        r_classes = set(rmat.rclass)
        if not supported_r_matrix_classes() & r_classes:
            if any(c.endswith('Matrix') for c in r_classes):
                # TODO(flying-sheep): #111 set stacklevel
                # https://github.com/theislab/anndata2ri/issues/111
                warn(f'Encountered Matrix class that is not supported: {r_classes}', stacklevel=2)
            return rmat
        for storage, mat_cls, idx, nnz in OPTIONS:
            if not supported_r_matrix_classes(storage=storage) & r_classes:
                continue
            coord_spec = idx(slots)
            data = (
                np.repeat(a=True, repeats=nnz(coord_spec))
                # we have pattern matrix without data (but always i and j!)
                if supported_r_matrix_classes(types='n') & r_classes
                else slots['x']
            )
            dtype = np.bool_ if supported_r_matrix_classes(types=('n', 'l')) & r_classes else np.float64
            return mat_cls((data, *coord_spec), shape=shape, dtype=dtype)

        msg = 'Should have hit one of the branches'
        raise AssertionError(msg)
