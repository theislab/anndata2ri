from functools import partial
from typing import Tuple, Type, Callable

import numpy as np
import pytest
from rpy2.rinterface import Sexp
from rpy2.robjects import baseenv, r, numpy2ri
from scipy import sparse

from anndata2ri import scipy2ri
from anndata2ri.rpy2_ext import importr
from anndata2ri.test_utils import conversions_rpy2py, ConversionModule


matrix = importr("Matrix")
methods = importr("methods")


dgc_empty = partial(methods.new, "dgCMatrix")
csc_empty = []

dgc = partial(r, "Matrix::Matrix((1:6)-3, 3, sparse = TRUE)")
csc_f = [[-2.0, 1.0], [-1.0, 2.0], [0.0, 3.0]]
dgr = partial(r, "as(matrix((1:6)-3, 2), 'dgRMatrix')")
csr_f = [[-2.0, 0.0, 2.0], [-1.0, 1.0, 3.0]]
dgt = partial(r, "Matrix::sparseMatrix(1:2, 3:2, x = 1:2, giveCsparse = FALSE)")
coo_f = [[0.0, 0.0, 1.0], [0.0, 2.0, 0.0]]

lgc = partial(r, "Matrix::Matrix(c(T, T, F, T, F, T), 2, sparse = TRUE)")
csc_b1 = [[1, 0, 0], [1, 1, 1]]
lgr = partial(r, "new('lgRMatrix', j = 1:0, p = c(0L, 0L, 1L, 2L), x = c(T, F), Dim = c(3L, 3L))")
csr_b1 = [[0, 0, 0], [0, 1, 0], [0, 0, 0]]
# lgt = ...
# coo_b1 = ...

ngc = partial(r, "as(Matrix::Matrix(c(T, T, F, T, F, T), 2, sparse = TRUE), 'nMatrix')")
csc_b2 = [[1, 0, 0], [1, 1, 1]]
# ngr = ...
# csr_b2 = ...
ngt = partial(r, "Matrix::sparseMatrix(1:2, 3:2, dims = c(3, 3), giveCsparse = FALSE)")
coo_b2 = [[0, 0, 1], [0, 1, 0], [0, 0, 0]]

mats = [
    pytest.param((0, 0), sparse.csc_matrix, np.floating, csc_empty, dgc_empty, id="dgC_empty"),
    pytest.param((3, 2), sparse.csc_matrix, np.floating, csc_f, dgc, id="dgC"),
    pytest.param((2, 3), sparse.csr_matrix, np.floating, csr_f, dgr, id="dgR"),
    pytest.param((2, 3), sparse.coo_matrix, np.floating, coo_f, dgt, id="dgT"),
    pytest.param((2, 3), sparse.csc_matrix, np.bool_, csc_b1, lgc, id="lgC"),
    pytest.param((3, 3), sparse.csr_matrix, np.bool_, csr_b1, lgr, id="lgR"),
    # pytest.param((?, ?), sparse.coo_matrix, np.bool_, coo_b1, lgt, id="lgT"),
    pytest.param((2, 3), sparse.csc_matrix, np.bool_, csc_b2, ngc, id="ngC"),
    # pytest.param((?, ?), sparse.csr_matrix, np.bool_, csr_b2, ngr, id="ngR"),
    pytest.param((3, 3), sparse.coo_matrix, np.bool_, coo_b2, ngt, id="ngT"),
]


@pytest.mark.parametrize("conversion", conversions_rpy2py)
@pytest.mark.parametrize("shape,cls,dtype,arr,dataset", mats)
def test_py2rpy(
    conversion: Callable[[ConversionModule, Callable[[], Sexp]], sparse.spmatrix],
    shape: Tuple[int, int],
    cls: Type[sparse.spmatrix],
    dtype: np.dtype,
    arr: np.ndarray,
    dataset: Callable[[], Sexp],
):
    sm = conversion(scipy2ri, dataset)
    assert isinstance(sm, cls)
    assert sm.shape == shape
    assert np.allclose(sm.toarray(), np.array(arr))

    dm = numpy2ri.converter.rpy2py(baseenv["as.matrix"](dataset()))
    assert np.allclose(sm.toarray(), dm)
