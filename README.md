# PGR-tk: A PanGenomic Research Took Kit

[![test_and_build](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml/badge.svg)](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml)

This repository is a project to provide Python and Rust libraries to facilitate pangenomics analysis. Several algorithms and data structures used for the Peregrine Genome Assembler are useful for Pangenomics analysis. This repo takes those algorithms and data structure, combining other handy 3rd party tools to expose them as a library in Python (with Rust code for those computing parts that need performance.) 

## Build

See `docker/Dockerfile.build_env-20.04` for a build enviroment under ubuntu 20.04.
With the proper build environment, just run `bash build.sh` to build all.

