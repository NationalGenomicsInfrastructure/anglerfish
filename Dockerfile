FROM continuumio/miniconda3

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish" \
      maintainer="remi-andre.olsen@scilifelab.se"

RUN apt-get update --allow-releaseinfo-change && apt-get install -y build-essential bash
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/anglerfish/bin:$PATH
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Add source files to the container
ADD . /usr/src/anglerfish
WORKDIR /usr/src/anglerfish

RUN pip install .
