import pytest
import numpy as np
from rpy2.robjects import baseenv, numpy2ri
from scipy import sparse

from anndata2ri import scipy2ri
from anndata2ri.test_utils import conversions_py2rpy

mats = [((2, 3), sparse.csr_matrix([[2.0, 0.0, 1.0], [0.0, 0.1, 0.0]]))]


@pytest.mark.parametrize("conversion", conversions_py2rpy)
@pytest.mark.parametrize("shape,dataset", mats)
def test_py2rpy_manual(conversion, shape, dataset):
    sm = conversion(scipy2ri, dataset)
    assert tuple(baseenv["dim"](sm)) == shape

    dm = numpy2ri.converter.py2rpy(baseenv["as.matrix"](sm))
    assert np.allclose(dataset.toarray(), dm)
