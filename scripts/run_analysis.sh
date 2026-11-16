#!/usr/bin/env bash
set -e   # stop at the first error

COUNTS=results/counts/count_matrix.tsv
META=metadata/samples.tsv
GTF=data/reference/annotation.gtf

Rscript scripts/01_sample_qc.R $COUNTS $META results/qc_test
Rscript scripts/02_DESeq2.R $COUNTS $META $GTF results/deseq2_test 0.05 1