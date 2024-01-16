FROM mambaorg/micromamba

LABEL author="Remi-Andre Olsen" \
      description="Anglerfish" \
      maintainer="remi-andre.olsen@scilifelab.se"

USER root

# Check for arm64 architecture to install minimap2
ARG TARGETARCH
RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Install compliation tools for minimap2
      apt-get update;\
      apt-get install -y curl build-essential libz-dev;\
    fi

ENV MINIMAP_VERSION=2.26

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Download minimap2
      curl -L https://github.com/lh3/minimap2/archive/refs/tags/v${MINIMAP_VERSION}.tar.gz | tar -zxvf - ;\
    fi

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      # Compile minimap2
      cd minimap2-${MINIMAP_VERSION};\
      make arm_neon=1 aarch64=1; \
      # Add to path
      mv minimap2 /usr/local/bin/;\
    fi

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /

RUN if [ "$TARGETARCH" = "arm64" ]; then \
      grep -v 'minimap2' /environment.yml > /environment.tmp.yml ;\
    else \
      cp /environment.yml /environment.tmp.yml ;\
    fi ;\
    chown $MAMBA_USER:$MAMBA_USER /environment.tmp.yml

# Add source files to the container
ADD . /usr/src/anglerfish
WORKDIR /usr/src/anglerfish

RUN micromamba install -y -n base -f /environment.tmp.yml && micromamba clean --all --yes

RUN eval "$(micromamba shell hook --shell bash)" && micromamba activate base && python -m pip install . 
USER $MAMBA_USER