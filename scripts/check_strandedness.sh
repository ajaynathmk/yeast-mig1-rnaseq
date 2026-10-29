#!/usr/bin/env bash
# Decide featureCounts -s (0/1/2) by counting how many reads are assigned to genes
# under each setting. Run on ONE sorted BAM after the first STAR alignment.
#
# Usage: scripts/check_strandedness.sh <sorted.bam> <annotation.gtf> [single|paired]
set -euo pipefail

BAM=${1:?bam required}
GTF=${2:?gtf required}
LAYOUT=${3:-single}

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

EXTRA=""
[ "$LAYOUT" = "paired" ] && EXTRA="-p --countReadPairs"

for s in 0 1 2; do
  featureCounts -T 4 -s "$s" $EXTRA -t exon -g gene_id -a "$GTF" -o "$TMP/s$s.txt" "$BAM" >/dev/null 2>&1
  awk -v s="$s" 'NR>1 { tot+=$2; if ($1=="Assigned") a=$2 }
                 END { printf "-s %d : %d / %d reads assigned (%.1f%%)\n", s, a, tot, 100*a/tot }' \
      "$TMP/s$s.txt.summary"
done

echo
echo "Reading the result:"
echo "  one of -s 1 / -s 2 is clearly higher (e.g. 85% vs 5%)  -> use that value"
echo "  -s 1 and -s 2 both about half of -s 0                  -> library is unstranded, use 0"
