
## Introduction

Thie directory contain docker file to build a docker image including the
pgr-tk and a jupyter lab server that can run some example notebook for
pangneome analysis.

This tutorual assumes the user is familiar with a typical Linux environment and docker.

## Build the Docker Image

After bullding the python wheel for version `0.x.y` (check the `../target/wheels` 
directory for proper version `0.x.y`).
The image uses python3.8.
 
```
cp ../target/wheels/pgrtk-0.x.y-cp38-cp38-linux_x86_64.whl .
docker build -t pgr-tk-ws:v0.x.y .
```

You can also use us a prebuilt docker image bye

```
docker pull cschin/pgr-tk-ws:v0.x.y
```
(check `https://hub.docker.com/r/cschin/pgr-tk-ws` the latest version.)


## Set up environment 

In a directory that you have write permission,

```
mkdir -p workdir
cd workdir
wget https://giab-data.s3.amazonaws.com/PGR-TK-Files/pgr-tk-HGRP-y1-evaluation-set-v0.tar
wget https://giab-data.s3.amazonaws.com/PGR-TK-Files/pgr-tk-example-code.zip
```

Untar the data file
```
tar xvf pgr-tk-HGRP-y1-evaluation-set-v0.tar
```

The tar ball contains the following data file  

```
data/
data/pgr-tk-HGRP-y1-evaluation-set-v0.agc # HPRC year 1 47 genomes (94 haplotype) + hg38 + hg19 + chm13 sequences in AGC format
data/pgr-tk-HGRP-y1-evaluation-set-v0.mdb # the SHIMMER index into the sequences 
data/pgr-tk-HGRP-y1-evaluation-set-v0.midx # auxilary index file for sequence names
data/pgr-tk-HGRP-y1-evaluation-set-v0_input # file used to generate the index 
data/AMY1A_gfa_view.png # AMY1A GFA example
```

Unzip the example notebooks

```
mkdir -p code && pushd code
unzip ../pgr-tk-example-code.zip
popd
```

Execute the Jupyter Lab server through docker

```
docker run -v $PWD:/wd/ -p 8888:8888 pgr-tk-ws:v0.x.y
```

or use a pre-built docker

Then follow the instruction from the Jupyter Lab output to connect to
the server from a browser.

For analyzing the whole 97 haplotyp human assembly, it is suggested
to have at least 64G RAM. You may use a remote server with enouge memory
connect to the server directly or through ssh tunneling.






