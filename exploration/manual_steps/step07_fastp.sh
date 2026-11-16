#!/usr/bin/env bash
# Step 7: adapter trimming with fastp (auto adapter detection for paired-end, minimum length 36).
set -euo pipefail
THREADS=${1:-4}
mkdir -p data/trimmed results/fastp
grep -q "data/trimmed/" .gitignore || echo "data/trimmed/" >> .gitignore
for r in $(cut -f1 metadata/samples.tsv | tail -n +2); do
  fastp -i data/raw/${r}_1.fastq.gz -I data/raw/${r}_2.fastq.gz \
        -o data/trimmed/${r}_1.fastq.gz -O data/trimmed/${r}_2.fastq.gz \
        --detect_adapter_for_pe --length_required 36 --thread "$THREADS" \
        --html results/fastp/${r}.html --json results/fastp/${r}.json \
        > results/fastp/${r}.log 2>&1
done
grep -E "reads passed filter|reads with adapter trimmed|Insert size peak" results/fastp/*.log
