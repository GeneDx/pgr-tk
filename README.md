# PGR-tk: A PanGenomic Research Took Kit

[![test_and_build](https://github.com/genedx/pgr-tk/actions/workflows/test_and_build.yml/badge.svg)](https://github.com/cschin/genedx/actions/workflows/test_and_build.yml)

This repository is a project to provide Python and Rust libraries to facilitate pangenomics analysis. Several algorithms and data structures used for the Peregrine Genome Assembler are useful for Pangenomics analysis as well. This repo takes those algorithms and data structure, combining other handy 3rd party tools to expose them as a library in Python (with Rust code for those computing parts that need performance.) 

## What is PGR-tk?

Research Preprint: 

[Multiscale Analysis of Pangenome Enables Improved Representation of Genomic Diversity For Repetitive And Clinically Relevant Genes](https://www.biorxiv.org/content/10.1101/2022.08.05.502980v2)

PGR-TK provides pangenome assembly management, query and Minimizer Anchored Pangenome (MAP) Graph Generation

![Pangenome Data Management and Minimizer Anchored Pangenome Graph Generation](/images/PGR_TK_Sketch_MAPG_construction.png)

With the MAP graph, we can use the "principal bundle decomposition" to study complicated structure variants and genome re-arragenment in the human populations.

![AMY1A Example](/images/AMY1A_example.png)


## Documentation, Usage and Examples

Command Line Tools:

PGR-TK provides the following tool to 

- create the PGR-TK sequence and index database
	-  `pgr-mdb`: create pgr minimizer database with AGC backend
	-  `pgr-make-frgdb`: create PGR-TK fragment minimizer database with frg format backend
- query the database to fetch sequences
	- `pgr-query`: query a PGR-TK pangenome sequence database, ouput the hit summary and generate fasta files from the target sequences
- generate MAP-graph in GFA format and principal bundle decomposition bed file
	- `pgr-pbundle-decomp`: generat the principal bundle decomposition though MAP Graph from a fasta file
- generate SVG from the principal bundle decomposition bed file
	- `pgr-pbundle-bed2svg`: generate SVG from a principal bundle bed file
- auxiliary tools
	- `pgr-pbundle-bed2sorted`: generate annotation file with a sorting order from the principal bundle decomposition
	- `pgr-pbundle-bed2dist`: generate alignment scores between sequences using bundle decomposition from a principal bundle bed file

For each comannd, `command --help` provides the detail usage information. 

The API documentation is at https://sema4-research.github.io/pgr-tk/

A collection of Jupyter Notebooks are at https://github.com/sema4-Research/pgr-tk-notebooks/

## Built Binaries

Check https://github.com/genedx/pgr-tk/releases


## Build

See `docker/Dockerfile.build_env-20.04` for a build enviroment under ubuntu 20.04.
With the proper build environment, just run `bash build.sh` to build all.

For example, on a Mac OS with Docker install, you can clone the repository and build a linux binary
within an Ubuntu 20.04 Linux distribution as follow:

1. Build the Docker image for a build environment:

```
git clone --recursive git@github.com:cschin/pgr-tk.git # clone the repo
cd pgr-tk/docker
ln -s Dockerfile.build_env-20.04 Dockerfile
docker build -t pgr-tk-build .
```

2. In the root directory of the repo `pgr-tk`:

Execute 
```
docker run -it --rm -v $PWD:/wd/pgr-tk pgr-tk-build /bin/bash 
```

3. Build the `pgr-tk` inside the docker container from the image `pgr-tk-build`

```
cd /wd/pgr-tk
bash build.sh
```

The build python wheels will be in `target/wheels` which can be installed for ubuntun 20.04 python3.8 distribution. You can install it in the `pgr-tk-build` image as well to test it out.

## Install stable verison v0.3.6 with Bioconda

If you have a conda install, you can try this to build an conda environment to use pgr-tk v0.3.6 (on linux only):

```
conda create -n pgr-tk python=3.8
conda activate pgr-tk
conda install -c bioconda -c conda-forge python_abi libstdcxx-ng=12 libclang13 pgr-tk=0.3.6
```

