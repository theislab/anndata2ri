|pypi| |conda| |rtd| |tests| |doi|

.. |pypi| image:: https://img.shields.io/pypi/v/anndata2ri
   :target: https://pypi.org/project/anndata2ri/
   :alt: PyPI Version

.. |conda| image:: https://img.shields.io/conda/vn/bioconda/anndata2ri
   :target: https://anaconda.org/bioconda/anndata2ri
   :alt: Bioconda Version

.. |rtd| image:: https://readthedocs.com/projects/icb-anndata2ri/badge/?version=latest&token=ee358f7efe36cbbd7d04db1b708fa81cefc44634ae7f3f8e0afcd03a1f0b1158
   :target: docs_
   :alt: Documentation Status

.. |tests| image:: https://github.com/theislab/anndata2ri/actions/workflows/run_tests.yml/badge.svg
   :target: https://github.com/theislab/anndata2ri/actions/workflows/run_tests.yml
   :alt: Unit Tests

.. |doi| image:: https://zenodo.org/badge/171714778.svg
   :target: https://zenodo.org/badge/latestdoi/171714778
   :alt: Publication DOI

AnnData ↭ SingleCellExperiment
==============================

RPy2 converter from AnnData_ to SingleCellExperiment_ and back. (For **details about conversion** see the docs_)

You can for example use it to process your data using both Scanpy_ and Seurat_, as described in this `example notebook`_

.. _AnnData: https://anndata.readthedocs.io/en/latest/
.. _SingleCellExperiment: http://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html
.. _docs: https://icb-anndata2ri.readthedocs-hosted.com/en/latest/
.. _Scanpy: https://scanpy.readthedocs.io/en/stable/
.. _Seurat: https://satijalab.org/seurat/
.. _`example notebook`: https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb

Installation
------------

.. code-block:: bash

   pip install anndata2ri
   # or
   conda install -c bioconda anndata2ri

Troubleshooting
---------------

If you have problems installing or importing anndata2ri,
please make sure you first:

1. Check the stack trace:
   If the error happens while installing or importing a dependency such as rpy2_,
   report your problem in that project’s bug tracker.
2. Search the issues_.
   At the time of writing 17 of the 29 bugs (60%) are invalid or rpy2 bugs / install problems.
3. Make sure you have a compatible R version: rpy2 requires R ≥ 3.6.

.. _rpy2: https://github.com/rpy2/rpy2#readme
.. _issues: https://github.com/theislab/anndata2ri/issues

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

Usage from Jupyter
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
