#!/usr/bin/env bash
# Step 8: FastQC on the trimmed reads.
set -euo pipefail
mkdir -p results/qc_trimmed
fastqc -t "${1:-4}" --outdir results/qc_trimmed data/trimmed/*.fastq.gz
