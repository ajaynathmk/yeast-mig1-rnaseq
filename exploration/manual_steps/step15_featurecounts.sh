#!/usr/bin/env bash
# Step 15: count read pairs per gene for all samples in one featureCounts run.
# Usage: scripts/step15_featurecounts.sh [threads] [strandedness, default 0]
set -euo pipefail
mkdir -p results/counts logs
BAMS=$(for r in $(cut -f1 metadata/samples.tsv | tail -n +2); do echo results/star/$r/Aligned.sortedByCoord.out.bam; done)
featureCounts -T "${1:-4}" -p --countReadPairs -s "${2:-0}" -t exon -g gene_id \
  -a data/reference/annotation.gtf -o results/counts/counts.txt $BAMS > logs/featurecounts.log 2>&1
