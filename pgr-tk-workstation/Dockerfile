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
        libomp5 \
    && rm -rf /var/lib/apt/lists/*


RUN . /opt/conda/bin/activate && \
conda install -y python=3.10 jupyterlab numpy networkx==2.4 matplotlib bokeh && conda clean -ya

#RUN conda install -y  --channel=conda-forge \
#        matplotlib \
#        numpy \
#    && conda clean -ya

RUN pip3 install -U --no-cache-dir \
    networkx==2.8.2 \
    papermill==2.3.4 \
    pydot==1.4.2 \
    scikit-learn

RUN echo deb http://http.us.debian.org/debian/ testing non-free contrib main > /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y libc6 libstdc++6
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /opt/conda/bin/../lib/libstdc++.so.6

COPY pgrtk-0.4.0-cp310-cp310-linux_x86_64.whl /tmp
RUN pip install /tmp/pgrtk-0.4.0-cp310-cp310-linux_x86_64.whl
RUN rm /tmp/pgrtk-0.4.0-cp310-cp310-linux_x86_64.whl
RUN mkdir -p /opt/bin/ /wd/
COPY jupyterlab.sh /opt/bin/
CMD /bin/bash /opt/bin/jupyterlab.sh
