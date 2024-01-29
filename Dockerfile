FROM mambaorg/micromamba as base

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish development version" \
      maintainer="remi-andre.olsen@scilifelab.se"

USER root

# Check for arm64 architecture to install minimap2
ARG TARGETARCH
ENV MINIMAP_VERSION=2.26
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Install compliation tools for minimap2
      apt-get update;\
      apt-get install -y curl build-essential libz-dev;\
      # Download minimap2
      curl -L https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP_VERSION}.tar.gz | tar -zxvf - ;\
      # Compile minimap2
      cd minimap2-${MINIMAP_VERSION};\
      make arm_neon=1 aarch64=1; \
      # Add to path
      mv minimap2 /usr/local/bin/;\
    fi

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /

# Remove minimap2 from environment.yml for arm64
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      grep -v 'minimap2' /environment.yml > /environment.tmp.yml ;\
    else \
      cp /environment.yml /environment.tmp.yml ;\
    fi ;\
    chown $MAMBA_USER:$MAMBA_USER /environment.tmp.yml

# Add source files to the container
ADD . /usr/src/anglerfish
WORKDIR /usr/src/anglerfish

# Activate the environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -n base -f /environment.tmp.yml && micromamba clean --all --yes

#####
# Devcontainer
#####
FROM base as devcontainer

# Useful tools for devcontainer
RUN apt-get update;\
    apt-get install -y git vim
RUN eval "$(micromamba shell hook --shell bash)" && python -m pip install -e .[dev]

#####
# Main
#####
FROM base as main
USER $MAMBA_USER
RUN eval "$(micromamba shell hook --shell bash)" && python -m pip install .[dev]