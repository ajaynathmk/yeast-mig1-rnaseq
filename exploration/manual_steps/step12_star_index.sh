#!/usr/bin/env bash
# Step 12: build the STAR index. sjdbOverhang = read length (151) - 1; SAindexNbases 10 suits the 12 Mb yeast genome.
set -euo pipefail
mkdir -p data/reference/star_index logs
STAR --runMode genomeGenerate \
     --runThreadN "${1:-4}" \
     --genomeDir data/reference/star_index \
     --genomeFastaFiles data/reference/genome.fa \
     --sjdbGTFfile data/reference/annotation.gtf \
     --sjdbOverhang 150 \
     --genomeSAindexNbases 10 \
     --outFileNamePrefix data/reference/star_index_ > logs/star_index.log 2>&1
