|pypi| |conda| |travis|

.. |pypi| image:: https://img.shields.io/pypi/v/anndata2ri
   :target: https://pypi.org/project/anndata2ri/

.. |conda| image:: https://img.shields.io/conda/vn/bioconda/anndata2ri
   :target: https://anaconda.org/bioconda/anndata2ri

.. |travis| image:: https://travis-ci.org/flying-sheep/anndata2ri.svg?branch=master
   :target: https://travis-ci.org/flying-sheep/anndata2ri

AnnData ↭ SingleCellExperiment
==============================

RPy2 converter from AnnData_ to SingleCellExperiment_ and back.

You can for example use it to process your data using both Scanpy_ and Seurat_, as described in this `example notebook`_

.. _AnnData: https://anndata.readthedocs.io/en/latest/
.. _SingleCellExperiment: http://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
.. _Scanpy: https://scanpy.readthedocs.io/en/stable/
.. _Seurat: https://satijalab.org/seurat/
.. _`example notebook`: https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb

Installation
------------

Install it directly from GitHub via ``pip`` (version 19.0 or higher).

.. code-block:: bash

   pip install git+https://github.com/flying-sheep/anndata2ri.git

You can install a locally checked out version with ``pip`` or ``flit``:

.. code-block:: bash

   cd anndata2ri
   pip install .
   # or
   flit install

Usage from Python
-----------------

Either use the converter manually …

.. code-block:: python

   import anndata2ri
   from rpy2.robjects import r
   from rpy2.robjects.conversion import localconverter

   with localconverter(anndata2ri.converter):
       adata = r('as(some_data, "SingleCellExperiment")')

… or activate it globally:

.. code-block:: python

   import anndata2ri
   from rpy2.robjects import r
   anndata2ri.activate()

   adata = r('as(some_data, "SingleCellExperiment")')

Usage from IPython
------------------
Activate the conversion before you load the extension:

.. code-block:: python

   import anndata2ri
   anndata2ri.activate()
   %load_ext rpy2.ipython

Now you can move objects from Python to R …

.. code-block:: python

   import scanpy.datasets as scd
   adata_paul = scd.paul15()

.. code-block:: r

   %%R -i adata_paul
   adata_paul  # class: SingleCellExperiment ...

… and back:

.. code-block:: r

   %%R -o adata_allen
   data(allen, package = 'scRNAseq')
   adata_allen <- as(allen, 'SingleCellExperiment')

.. code-block:: python

   print(adata_allen)  # AnnData object with ...
