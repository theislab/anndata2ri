from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pytest
from rpy2.robjects import baseenv, numpy2ri, r
from scipy import sparse

from anndata2ri import scipy2ri
from anndata2ri._rpy2_ext import importr


if TYPE_CHECKING:
    from collections.abc import Callable

    from rpy2.rinterface import Sexp

    from anndata2ri.test_utils import R2Py


matrix = importr('Matrix')
methods = importr('methods')


dgc_empty = partial(methods.new, 'dgCMatrix')
csc_empty = []

dgc = partial(r, 'Matrix::Matrix((1:6)-3, 3, sparse = TRUE)')
csc_f = [[-2.0, 1.0], [-1.0, 2.0], [0.0, 3.0]]
dgr = partial(r, "as(matrix((1:6)-3, 2), 'dgRMatrix')")
csr_f = [[-2.0, 0.0, 2.0], [-1.0, 1.0, 3.0]]
dgt = partial(r, 'Matrix::sparseMatrix(1:2, 3:2, x = 1:2, giveCsparse = FALSE)')
coo_f = [[0.0, 0.0, 1.0], [0.0, 2.0, 0.0]]

lgc = partial(r, 'Matrix::Matrix(c(T, T, F, T, F, T), 2, sparse = TRUE)')
csc_b1 = [[1, 0, 0], [1, 1, 1]]
lgr = partial(r, "new('lgRMatrix', j = 1:0, p = c(0L, 0L, 1L, 2L), x = c(T, F), Dim = c(3L, 3L))")
csr_b1 = [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
lgt = partial(r, "Matrix::sparseMatrix(1:2, 3:2, x = c(TRUE, FALSE), repr = 'T')")
coo_b1 = [[0, 0, 1], [0, 0, 0]]

ngc = partial(r, "as(Matrix::Matrix(c(T, T, F, T, F, T), 2, sparse = TRUE), 'nMatrix')")
csc_b2 = [[1, 0, 0], [1, 1, 1]]
ngr = partial(r, "Matrix::sparseMatrix(1:2, 3:2, dims = c(3, 3), repr = 'R')")
csr_b2 = [[0, 0, 1], [0, 1, 0], [0, 0, 0]]
ngt = partial(r, 'Matrix::sparseMatrix(1:2, 3:2, dims = c(3, 3), giveCsparse = FALSE)')
coo_b2 = [[0, 0, 1], [0, 1, 0], [0, 0, 0]]

ddi = partial(r, 'Matrix::Diagonal(x = c(1, -2, 3))')
dia_f = [[1.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 0.0, 3.0]]
ddi_unit = partial(r, 'Matrix::Diagonal(3)')
dia_f_unit = np.eye(3)
ldi = partial(r, 'Matrix::Diagonal(x = c(TRUE, FALSE, TRUE))')
dia_b = [[1, 0, 0], [0, 0, 0], [0, 0, 1]]

mats = [
    pytest.param((0, 0), sparse.csc_array, np.float64, csc_empty, dgc_empty, id='dgC_empty'),
    pytest.param((3, 2), sparse.csc_array, np.float64, csc_f, dgc, id='dgC'),
    pytest.param((2, 3), sparse.csr_array, np.float64, csr_f, dgr, id='dgR'),
    pytest.param((2, 3), sparse.coo_array, np.float64, coo_f, dgt, id='dgT'),
    pytest.param((2, 3), sparse.csc_array, np.bool_, csc_b1, lgc, id='lgC'),
    pytest.param((3, 3), sparse.csr_array, np.bool_, csr_b1, lgr, id='lgR'),
    pytest.param((2, 3), sparse.coo_array, np.bool_, coo_b1, lgt, id='lgT'),
    pytest.param((2, 3), sparse.csc_array, np.bool_, csc_b2, ngc, id='ngC'),
    pytest.param((3, 3), sparse.csr_array, np.bool_, csr_b2, ngr, id='ngR'),
    pytest.param((3, 3), sparse.coo_array, np.bool_, coo_b2, ngt, id='ngT'),
    pytest.param((3, 3), sparse.dia_array, np.float64, dia_f, ddi, id='ddi'),
    pytest.param((3, 3), sparse.dia_array, np.float64, dia_f_unit, ddi_unit, id='ddi_unit'),
    pytest.param((3, 3), sparse.dia_array, np.bool_, dia_b, ldi, id='ldi'),
]


@pytest.mark.parametrize(('shape', 'cls', 'dtype', 'arr', 'dataset'), mats)
def test_mats(
    *,
    r2py: R2Py,
    shape: tuple[int, int],
    cls: type[sparse.sparray],
    dtype: np.dtype,
    arr: np.ndarray,
    dataset: Callable[[], Sexp],
) -> None:
    sm = r2py(scipy2ri, dataset)
    assert isinstance(sm, cls)
    assert sm.shape == shape
    assert sm.dtype == dtype
    assert np.allclose(sm.toarray(), np.array(arr))

    dm = numpy2ri.converter.rpy2py(baseenv['as.matrix'](dataset()))
    assert np.allclose(sm.toarray(), dm)
