.. _home-page-install:

**************
How to install
**************

.. autosummary::
   :toctree: generated

Please install nextflow `Nextflow <https://www.nextflow.io/>`_ and `Singularity <https://sylabs.io/>`_ or `Docker <https://www.docker.com/>`_ before.

Then download the repo:

.. code-block:: console

  git clone --depth 1 --recurse-submodules https://github.com/biocorecrg/MOP2.git


You can use **INSTALL.sh** to download the **guppy 3.4.5** or you can download the version you prefer by adding an extra parameter. 

.. note::
  
  Please consider that the support of VBZ compression of fast5 started with version 3.4.X. 


.. code-block:: console
  
  cd MoP2; bash INSTALL.sh 

or

.. code-block:: console

  cd MoP2; bash INSTALL.sh 4.0.15
 
 
Testing
============

.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -with-singularity -bg -profile local > log

.. tip::

  You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.


