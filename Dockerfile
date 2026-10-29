FROM mambaorg/micromamba:1.5.8

# Install everything from environment.yml into the base env of the image
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Makes `docker build` RUN steps (and later commands) run inside the env
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /work
CMD ["/bin/bash"]
