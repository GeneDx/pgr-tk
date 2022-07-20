# PGR-tk: A PanGenomic Research Took Kit

[![test_and_build](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml/badge.svg)](https://github.com/cschin/pgr-tk/actions/workflows/test_and_build.yml)

This repository is a project to provide Python and Rust libraries to facilitate pangenomics analysis. Several algorithms and data structures used for the Peregrine Genome Assembler are useful for Pangenomics analysis. This repo takes those algorithms and data structure, combining other handy 3rd party tools to expose them as a library in Python (with Rust code for those computing parts that need performance.) 

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


