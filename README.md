# MoP2- DSL2 version of Master of Pores

![MOP2](https://github.com/biocorecrg/MoP2/blob/master/img/master_red.jpg?raw=true)


Inspired by Metallica's [Master Of Puppets](https://www.youtube.com/watch?v=S7blkui3nQc)

## Install
Please install nextflow and singularity or docker before.

Then download the repo:

```
git clone --depth 1 --recurse-submodules git@github.com:biocorecrg/MoP2.git
```

You can use INSTALL.sh to download the version 3.4.5 of guppy or you can replace it with the version you prefer. Please consider that the support of VBZ compression of fast5 started with version 3.4.X. 

```
cd MoP2; sh INSTALL.sh
```

## Testing
You can replace ```-with-singularity``` with ```-with-docker``` if you want to use the docker engine.

```
cd mop_preprocess
nextflow run mop_preprocess.nf -with-singularity -bg -profile local > log

```

## Documentation

In progress... [here](https://biocorecrg.github.io/MoP2/docs/)
