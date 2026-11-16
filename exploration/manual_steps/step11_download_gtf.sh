#!/usr/bin/env bash
# Step 11: download the gene annotation GTF (Ensembl Genomes fungi release 56).
set -euo pipefail
mkdir -p data/reference
wget -O data/reference/annotation.gtf.gz http://ftp.ensemblgenomes.org/pub/fungi/release-56/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.56.gtf.gz
gzip -t data/reference/annotation.gtf.gz && gunzip -f data/reference/annotation.gtf.gz
