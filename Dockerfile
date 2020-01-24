FROM continuumio/miniconda3

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish" \
      maintainer="remi-andre.olsen@scilifelab.se"

RUN apt-get update && apt-get install -y build-essential bash
RUN conda update -n base -c defaults conda && conda config --add channels conda-forge && conda config --add channels bioconda
RUN conda create -n anglerfish && conda update --all --yes
RUN conda install -n anglerfish --yes multiqc==1.8 bioconda::minimap2 bioconda::fastqc biopython python-levenshtein
ENV PATH /opt/conda/envs/anglerfish/bin:$PATH
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Add source files to the container
ADD . /usr/src/anglerfish
WORKDIR /usr/src/anglerfish

RUN pip install .
