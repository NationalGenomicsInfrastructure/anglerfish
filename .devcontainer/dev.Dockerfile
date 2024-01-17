FROM mambaorg/micromamba

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish development version" \
      maintainer="remi-andre.olsen@scilifelab.se"

USER root

# Useful tools for devcontainer
RUN apt-get update;\
    apt-get install -y git

# Check for arm64 architecture to install minimap2
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Install compliation tools for minimap2
      apt-get update;\
      apt-get install -y curl build-essential libz-dev;\
    fi


USER $MAMBA_USER
ENV MINIMAP_VERSION=2.26

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Download minimap2
      curl -L https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP_VERSION}.tar.gz | tar -zxvf - ;\
    fi

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Compile minimap2
      cd minimap2-${MINIMAP_VERSION};\
      make arm_neon=1 aarch64=1; \
      # Add to PATH
      echo "export PATH=$(pwd):\$PATH" >> ~/.bashrc ;\
    fi

# Add requirements files to the container
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /tmp/

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      grep -v 'minimap2' /tmp/environment.yml > /tmp/environment.tmp.yml ;\
    else \
      cp /tmp/environment.yml /tmp/environment.tmp.yml ;\
    fi

# Activate the environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN micromamba install -y -f /tmp/environment.tmp.yml && micromamba clean --all --yes