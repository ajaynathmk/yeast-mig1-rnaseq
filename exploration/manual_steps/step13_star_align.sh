#!/usr/bin/env bash
# Step 13: align trimmed reads with STAR -> coordinate-sorted BAM per sample.
# --alignIntronMax/--alignMatesGapMax 2000: yeast introns are short; --limitBAMsortRAM avoids the BAM-sorting memory error.
set -euo pipefail
THREADS=${1:-4}
mkdir -p logs
for r in $(cut -f1 metadata/samples.tsv | tail -n +2); do
  mkdir -p results/star/$r
  rm -rf results/star/$r/_STARtmp
  STAR --runThreadN "$THREADS" \
       --genomeDir data/reference/star_index \
       --readFilesIn data/trimmed/${r}_1.fastq.gz data/trimmed/${r}_2.fastq.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix results/star/$r/ \
       --outSAMtype BAM SortedByCoordinate \
       --alignIntronMax 2000 --alignMatesGapMax 2000 \
       --limitBAMsortRAM 3000000000 > logs/star_$r.log 2>&1
  grep -E "Uniquely mapped reads %|Number of input reads" results/star/$r/Log.final.out | sed "s/^/$r /"
done
