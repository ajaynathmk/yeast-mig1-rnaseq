#!/usr/bin/env bash
# Step 9: MultiQC report combining trimmed-read FastQC and fastp results.
set -euo pipefail
multiqc results/qc_trimmed results/fastp -o results/multiqc_trimmed -f
