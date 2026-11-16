#!/usr/bin/env bash
# Step 18: DESeq2 (KO vs WT), annotation, volcano plot and GO enrichment.
# Usage: scripts/step18_deseq2.sh [padj_cutoff=0.05] [lfc_cutoff=1]
set -euo pipefail
mkdir -p results/deseq2 logs
Rscript scripts/02_DESeq2.R results/counts/count_matrix.tsv metadata/samples.tsv \
  data/reference/annotation.gtf results/deseq2 "${1:-0.05}" "${2:-1}" 2>&1 | tee logs/deseq2.log
