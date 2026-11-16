#!/usr/bin/env bash
# Step 16: turn counts.txt into a clean gene x sample matrix (columns named by run accession).
set -euo pipefail
tail -n +2 results/counts/counts.txt | cut -f1,7- \
  | sed '1s#results/star/\(SRR[0-9]*\)/Aligned.sortedByCoord.out.bam#\1#g' \
  > results/counts/count_matrix.tsv
head -3 results/counts/count_matrix.tsv
