FROM mambaorg/micromamba

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish development version" \
      maintainer="remi-andre.olsen@scilifelab.se"

USER root
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /
RUN micromamba env create -n anglerfish && micromamba install -y -n anglerfish -f /environment.yml && micromamba clean --all --yes 
ENV PATH /opt/conda/envs/anglerfish/bin:$PATH

# Add source files to the container
ADD . /usr/src/anglerfish
WORKDIR /usr/src/anglerfish
RUN eval "$(micromamba shell hook --shell bash)" && micromamba activate anglerfish && python -m pip install .[dev]
USER $MAMBA_USER
