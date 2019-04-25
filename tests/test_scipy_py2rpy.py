import pytest
import numpy as np
from rpy2.robjects import baseenv, numpy2ri
from scipy import sparse

from anndata2ri import scipy2ri
from anndata2ri.test_utils import conversions_py2rpy


mats = [
    pytest.param((0, 0), sparse.csr_matrix((0, 0)), "gR", id="csr-empty"),
    pytest.param((2, 3), sparse.csr_matrix([[2.0, 0.0, 1.0], [0.0, 0.1, 0.0]]), "gR", id="csr"),
    pytest.param((3, 2), sparse.csc_matrix([[2.0, 0.0], [1.0, 0.0], [0.1, 0.0]]), "gC", id="csc"),
    pytest.param((2, 4), sparse.coo_matrix([[2.0, 0.0, 1.0, 0.0], [0.0, 0.1, 0.0, 3.0]]), "gT", id="coo"),
    pytest.param((4, 4), sparse.dia_matrix(([2.0, 0.4, 1.0, 0.0], [0]), (4, 4)), "di", id="dia"),
    pytest.param((3, 3), sparse.dia_matrix(([0.0, 0.0, 0.0], [0]), (3, 3)), "di", id="dia_empty"),
    pytest.param((4, 4), sparse.dia_matrix(([1.0, 1.0, 1.0, 1.0], [0]), (4, 4)), "di", id="dia_unit"),
]


@pytest.mark.parametrize("typ", ["l", "d"])
@pytest.mark.parametrize("conversion", conversions_py2rpy)
@pytest.mark.parametrize("shape,dataset,cls", mats)
def test_py2rpy(typ, conversion, shape, dataset, cls):
    if typ == "l":
        dataset = dataset.astype(bool)
    sm = conversion(scipy2ri, dataset)
    assert f"{typ}{cls}Matrix" in set(sm.rclass)
    assert tuple(baseenv["dim"](sm)) == shape

    dm = numpy2ri.converter.py2rpy(baseenv["as.matrix"](sm))
    assert np.allclose(dataset.toarray(), dm)
