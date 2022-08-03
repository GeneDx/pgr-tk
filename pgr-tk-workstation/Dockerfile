FROM continuumio/miniconda3:latest

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        curl \
        git \
        graphviz \
        graphviz-dev \
        gawk \
        minimap2 \
        samtools \
        time \
        wget \
        pdf2svg \
        awscli \
        vim \
        build-essential \
        zlib1g \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y  --channel=conda-forge \
        matplotlib \
        numpy \
    && conda clean -ya

RUN . /opt/conda/bin/activate && \
conda install -y python=3.8.10 jupyterlab numpy networkx==2.4 matplotlib bokeh && conda clean -ya

RUN pip3 install -U --no-cache-dir \
    networkx==2.8.2 \
    papermill==2.3.4 \
    pydot==1.4.2 \
    sklearn

COPY pgrtk-0.3.4-cp38-cp38-linux_x86_64.whl /tmp
RUN pip install /tmp/pgrtk-0.3.4-cp38-cp38-linux_x86_64.whl
RUN rm /tmp/pgrtk-0.3.4-cp38-cp38-linux_x86_64.whl
RUN mkdir -p /opt/bin/ /wd/
COPY jupyterlab.sh /opt/bin/
CMD /bin/bash /opt/bin/jupyterlab.sh
