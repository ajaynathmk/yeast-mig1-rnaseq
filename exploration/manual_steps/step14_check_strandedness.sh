#!/usr/bin/env bash
# Step 14: decide featureCounts -s (0/1/2) by comparing % of read pairs assigned under each setting.
# Unstranded libraries give ~80% for -s 0 and ~46% for -s 1 and -s 2 (what we saw for all four samples).
set -euo pipefail
mkdir -p results/strand_test
for r in $(cut -f1 metadata/samples.tsv | tail -n +2); do
  for s in 0 1 2; do
    featureCounts -T "${1:-4}" -p --countReadPairs -s $s -t exon -g gene_id \
        -a data/reference/annotation.gtf \
        -o results/strand_test/${r}_s$s.txt \
        results/star/$r/Aligned.sortedByCoord.out.bam > /dev/null 2>&1
    awk -v r=$r -v s=$s 'NR>1{t+=$2; if($1=="Assigned") a=$2} END{printf "%s  -s %d : %.1f%% of read pairs assigned to genes\n", r, s, 100*a/t}' results/strand_test/${r}_s$s.txt.summary
  done
done
