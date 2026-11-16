FROM mambaorg/micromamba:1.5.8-noble

# 1. Install the exact tool versions into the container's base environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.lock.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 2. Copy the workflow (code only, no data)
WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER Snakefile config.yaml ./
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ scripts/
COPY --chown=$MAMBA_USER:$MAMBA_USER metadata/ metadata/

# 3. Let any user write inside /pipeline (the container is run as the host user)
USER root
RUN chmod -R a+rwX /pipeline
USER $MAMBA_USER

# 4. Running the container runs the pipeline
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "snakemake"]
CMD ["--cores", "4"]