# Install

Please install nextflow and singularity or docker before.

Then download the repo:

.. code-block:: console

  git clone --depth 1 --recurse-submodules git@github.com:biocorecrg/MoP2.git


You can use **INSTALL.sh** to download the **guppy 3.4.5** or you can replace it with the version you prefer. 

.. note::
  
  Please consider that the support of VBZ compression of fast5 started with version 3.4.X. 


.. code-block:: console
  
  cd MoP2; sh INSTALL.sh


## Testing

.. code-block:: console

  cd mop_preprocess

  nextflow run mop_preprocess.nf -with-singularity -bg -profile local > log

.. tip::

  You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.



