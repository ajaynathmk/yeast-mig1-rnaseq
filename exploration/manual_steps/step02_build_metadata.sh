#!/usr/bin/env bash
# Step 2: fetch the ENA run table for PRJNA412405 and build the sample sheet (2 WT vs 2 mig1 KO),
# the FASTQ URL list and the expected FASTQ sizes. Run from the project root.
set -euo pipefail
mkdir -p metadata
ACC=PRJNA412405

curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=run_accession,experiment_accession,sample_alias,sample_title,library_layout,read_count,fastq_ftp&format=tsv" > metadata/ena_raw.tsv

{ printf 'run_accession\tsample\tcondition\treplicate\tfastq_ftp\n'
  awk -F'\t' 'NR>1 && ($4=="ScY01a"||$4=="mig1"){print $1"\t"$4"\t"($4=="mig1"?"KO":"WT")"\t"$7}' metadata/ena_raw.tsv \
  | sort -k3,3 -k1,1 \
  | awk -F'\t' 'BEGIN{OFS="\t"}{n[$3]++; print $1,$2,$3,n[$3],$4}'
} > metadata/samples.tsv

# one URL per line (ftp:// prefix; download_one.sh tries https, http and ftp anyway)
cut -f5 metadata/samples.tsv | tail -n +2 | tr ';' '\n' | sed 's#^#ftp://#' > metadata/fastq_urls.txt

# expected file sizes, used to decide that a download is complete
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes&format=tsv" > metadata/fastq_sizes_raw.tsv
awk -F'\t' 'NR>1{n=split($2,u,";"); split($3,b,";"); for(i=1;i<=n;i++){k=split(u[i],p,"/"); print p[k]"\t"b[i]}}' metadata/fastq_sizes_raw.tsv > metadata/fastq_sizes.tsv

column -t -s $'\t' metadata/samples.tsv
