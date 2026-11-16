#!/usr/bin/env bash
# Step 1: create the micromamba environment with every tool used in the workflow.
set -euo pipefail
micromamba create -y -n yeast-rnaseq -c conda-forge -c bioconda python=3.11 fastqc star subread samtools
micromamba install -y -n yeast-rnaseq -c conda-forge -c bioconda multiqc "snakemake=7.32.4" "pulp=2.7.0" pandas requests wget fastp
micromamba install -y -n yeast-rnaseq -c conda-forge -c bioconda "r-base>=4.3" r-ggplot2 r-ggrepel r-pheatmap bioconductor-deseq2 bioconductor-apeglm
micromamba install -y -n yeast-rnaseq -c conda-forge -c bioconda r-tidyverse r-rcolorbrewer r-rmarkdown bioconductor-rtracklayer bioconductor-enhancedvolcano
micromamba env export -n yeast-rnaseq --from-history > environment.yml
