|travis|

.. |travis| image:: https://travis-ci.org/flying-sheep/anndata2ri.svg?branch=master
   :target: https://travis-ci.org/flying-sheep/anndata2ri

AnnData ↭ SingleCellExperiment
==============================

RPy2 converter from AnnData to SCE and back.

Usage from Python
-----------------

Either use the converter manually …

.. code-block:: python

   import anndata2ri
   from rpy2.robjects import r
   from rpy2.robjects.conversion import ConversionContext

   with ConversionContext(anndata2ri.create_converter()):
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
