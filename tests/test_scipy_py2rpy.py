from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pytest
from rpy2.robjects import baseenv, numpy2ri
from scipy import sparse

from anndata2ri import scipy2ri


if TYPE_CHECKING:
    from anndata2ri.test_utils import Py2R


mats = [
    pytest.param((0, 0), sparse.csr_matrix((0, 0)), 'gR', id='csr-empty'),
    pytest.param((2, 3), sparse.csr_matrix([[2.0, 0.0, 1.0], [0.0, 0.1, 0.0]]), 'gR', id='csr'),
    pytest.param((3, 2), sparse.csc_matrix([[2.0, 0.0], [1.0, 0.0], [0.1, 0.0]]), 'gC', id='csc'),
    pytest.param((2, 4), sparse.coo_matrix([[2.0, 0.0, 1.0, 0.0], [0.0, 0.1, 0.0, 3.0]]), 'gT', id='coo'),
    pytest.param((4, 4), sparse.dia_matrix(([2.0, 0.4, 1.0, 0.0], [0]), (4, 4)), 'di', id='dia'),
    pytest.param((3, 3), sparse.dia_matrix(([0.0, 0.0, 0.0], [0]), (3, 3)), 'di', id='dia_empty'),
    pytest.param((4, 4), sparse.dia_matrix(([1.0, 1.0, 1.0, 1.0], [0]), (4, 4)), 'di', id='dia_unit'),
]


@pytest.mark.parametrize('typ', ['l', 'd'])
@pytest.mark.parametrize(('shape', 'dataset', 'cls'), mats)
def test_mats(
    py2r: Py2R,
    typ: Literal['l', 'd'],
    shape: tuple[int, ...],
    dataset: sparse.spmatrix,
    cls: str,
) -> None:
    if typ == 'l':
        dataset = dataset.astype(bool)
    sm = py2r(scipy2ri, dataset)
    assert f'{typ}{cls}Matrix' in set(sm.rclass)
    assert tuple(baseenv['dim'](sm)) == shape

    dm = numpy2ri.converter.py2rpy(baseenv['as.matrix'](sm))
    assert np.allclose(dataset.toarray(), dm)
