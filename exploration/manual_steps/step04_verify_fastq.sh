#!/usr/bin/env bash
# Step 4: check every FASTQ is a valid gzip and that R1/R2 read counts equal ENA's read_count.
set -uo pipefail
status=0
for r in $(cut -f1 metadata/samples.tsv | tail -n +2); do
  gzip -t data/raw/${r}_1.fastq.gz data/raw/${r}_2.fastq.gz || { echo "$r gzip FAILED"; status=1; continue; }
  exp=$(awk -F'\t' -v r=$r '$1==r{print $6}' metadata/ena_raw.tsv)
  n1=$(( $(zcat data/raw/${r}_1.fastq.gz | wc -l) / 4 ))
  n2=$(( $(zcat data/raw/${r}_2.fastq.gz | wc -l) / 4 ))
  if [ "$n1" = "$exp" ] && [ "$n2" = "$exp" ]; then echo "$r OK  expected=$exp R1=$n1 R2=$n2"
  else echo "$r MISMATCH expected=$exp R1=$n1 R2=$n2"; status=1; fi
done
exit $status
