#!/usr/bin/env bash
# Step 5: FastQC on the raw reads.
set -euo pipefail
mkdir -p results/qc
fastqc -t "${1:-4}" --outdir results/qc data/raw/*.fastq.gz
