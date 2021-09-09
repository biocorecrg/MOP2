## Version 2.0
* Completely rewritten using the powerful [DSL2](https://www.nextflow.io/docs/latest/dsl2.html).
* Subworkflows are stored in the independent repository [BioNextflow](https://github.com/biocorecrg/BioNextflow).
* mop_preprocess (formerly known as nanoPreprocess + nanoPreprocessSimple)
  * now can read multiple runs per time
  * can demultiplex fast5 using guppy too
  * deeplexicon can be run on GPU too
  * Parameters of each tool are stored in a tsv file. We have three different ones already pre-set for cDNA, DNA and dRNA (option --pars_tools)
  * demultiplexing, filtering, mapping and counting can be switched off by setting "NO" as a parameter
  * saveSpace can be set to "YES" to reduce the amount of disk space required. WARNING This will prevent the possibility to resume!
  * Merged old NanoPreprocess and NanoPreprocessSimple in a mop_preprocess. Using fastq or fast5 will switch among the two executions.

* mop_mod (formerly known as nanoMod)
  * now you can launch each analysis independently
  * 4 workflows: epinano, nanopolish + nanocompore, tombo model_sample_compare and tombo level_sample_compare  

* ...

## Version 1.1
* Added a new module called NanoPreprocessSimple that starts from fastq files instead of fast5 files. It allows the analysis of multiple files at a time.
* Added support to vbz compressed fast5 https://github.com/nanoporetech/vbz_compression in NanoPreprocess, NanoMod and NanoTail
* NanoPreprocess now outputs also CRAM files and can do downsampling with the parameter --downsampling
* NanoPreprocess allows performing variant calling using medaka (BETA)
* NanoPreprocess allows performing demultiplexing with GUPPY
* Added plots for Epinano output in NanoMod
* Added a conversion of Tombo results in bed format in NanoMod
* Added a INSTALL.sh file for automatically retrieve guppy 3.4.5 from https://mirror.oxfordnanoportal.com/, place it in NanoPreprocess/bin and making the required links
* Added profiles for being used locally and on the CRG SGE cluster


## Version 1.0
This is the original version published in the paper [MasterOfPores: A Workflow for the Analysis of Oxford Nanopore Direct RNA Sequencing Datasets](https://www.frontiersin.org/articles/10.3389/fgene.2020.00211/full)
