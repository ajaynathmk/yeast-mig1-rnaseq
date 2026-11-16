#!/usr/bin/env bash
# Step 10: download the S. cerevisiae R64-1-1 genome FASTA (Ensembl Genomes fungi release 56).
set -euo pipefail
mkdir -p data/reference
wget -O data/reference/genome.fa.gz http://ftp.ensemblgenomes.org/pub/fungi/release-56/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gzip -t data/reference/genome.fa.gz && gunzip -f data/reference/genome.fa.gz
