import numpy as np
from anndata import AnnData
from rpy2.robjects import conversion, baseenv

import anndata2ri


def test_py2rpy_simple():
    ad = AnnData(np.array([[1, 2, 3], [0.3, 0.2, 0.1]]), dict(cluster=[1, 2]))
    with conversion.ConversionContext(anndata2ri.create_converter()) as c:
        ex = c.py2rpy(ad)
        assert baseenv["dim"](ex).tolist() == [3, 2]
