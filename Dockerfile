FROM frolvlad/alpine-miniconda3

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish" \
      maintainer="remi-andre.olsen@scilifelab.se"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/anglerfish/bin:$PATH

# Add source files to the container
ADD . /usr/src/anglerfish
WORKDIR /usr/src/anglerfish

RUN python -m pip install .
