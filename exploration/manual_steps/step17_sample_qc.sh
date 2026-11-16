#!/usr/bin/env bash
# Step 17: sample-level QC (PCA, sample distances, correlations) before any gene testing.
set -euo pipefail
mkdir -p results/sample_qc logs
Rscript scripts/01_sample_qc.R results/counts/count_matrix.tsv metadata/samples.tsv results/sample_qc 2>&1 | tee logs/sample_qc.log
