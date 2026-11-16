#!/usr/bin/env bash
# Step 6: MultiQC report of the raw-read FastQC results.
set -euo pipefail
multiqc results/qc -o results/multiqc_raw -f
