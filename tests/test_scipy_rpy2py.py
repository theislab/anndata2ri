import numpy as np
import pytest
from rpy2.robjects import baseenv, r, numpy2ri
from rpy2.robjects.packages import importr
from scipy import sparse

from anndata2ri import scipy2ri
from anndata2ri.test_utils import conversions_rpy2py


matrix = importr("Matrix")
methods = importr("methods")


mats = [
    pytest.param((0, 0), lambda: methods.new("dgCMatrix"), sparse.csc_matrix, id="dsC_empty"),
    pytest.param((3, 2), lambda: r("Matrix::Matrix((1:6)-3, 3, sparse = TRUE)"), sparse.csc_matrix, id="dsC"),
    # TODO: moar
]


@pytest.mark.parametrize("conversion", conversions_rpy2py)
@pytest.mark.parametrize("shape,dataset,cls", mats)
def test_py2rpy(conversion, shape, dataset,cls):
    sm = conversion(scipy2ri, dataset)
    assert isinstance(sm, cls)
    assert sm.shape == shape

    dm = numpy2ri.converter.rpy2py(baseenv["as.matrix"](dataset()))
    assert np.allclose(sm.toarray(), dm)
