# Yeast MIG1 deletion RNA-seq pipeline

Reproducible Snakemake + Docker workflow for differential expression in *Saccharomyces cerevisiae*:
FASTQ download (ENA) → FastQC → STAR → featureCounts → MultiQC → DESeq2 → GO / GSEA / KEGG.
Compares a *mig1*Δ deletion strain with wild type. The gene is a config setting, so the same
workflow runs for other deletion strains in the same study.

## Data

- GEO series: GSE104356 (BioProject: `PRJNA…` ← fill in), wild type vs *mig1*Δ, grown at 40 °C
- Design: 2 vs 2 biological replicates (low power, see Limitations)
- Citation: add the dataset publication listed on the GEO page

## Quick start

```bash
git clone https://github.com/ajaynathmk/yeast-mig1-rnaseq-pipeline.git
cd yeast-mig1-rnaseq-pipeline

docker build -t yeast-rnaseq .
docker run -it --rm --user "$(id -u):$(id -g)" -v "$PWD":/work yeast-rnaseq bash

# inside the container
python scripts/build_metadata.py PRJNA… --ko-pattern mig1   # then check metadata/samples.tsv by hand
snakemake -n --cores 4                                       # dry run
snakemake --cores 4                                          # full run
```

Check library strandedness once (see `scripts/check_strandedness.sh`) and set `strandedness` in `config.yaml`.

## Workflow

| Stage | Tool | Output |
|---|---|---|
| QC | FastQC, MultiQC | `results/multiqc/multiqc_report.html` |
| Alignment | STAR (R64-1-1) | `results/star/<run>/Aligned.sortedByCoord.out.bam` |
| Counting | featureCounts | `results/counts/counts.txt` |
| Differential expression | DESeq2 (apeglm shrinkage) | `results/deseq2/` |
| Enrichment | clusterProfiler | `results/enrichment/` |

## Results

_Add your findings here: MultiQC summary, PCA, volcano plot, number of up/down genes, the MIG1 positive-control check,
top enriched GO terms._

## Limitations

- n = 2 per group; modest effects may not reach significance
- Comparison is made under heat stress, not standard glucose-repression conditions
- Enrichment results are hypothesis-generating

## Repository layout

```
config.yaml           # gene, strandedness, thresholds, reference URLs
Snakefile             # workflow
scripts/              # metadata builder, strandedness check, DESeq2, enrichment
metadata/samples.tsv  # sample sheet (verified by hand)
results/              # figures and tables
```

## Credits and citations

DESeq2, STAR, featureCounts, apeglm, clusterProfiler, Snakemake, and the dataset publication.
Project structure inspired by the public repository *Dockerized_Yeast_Rnaseq_pipeline*.
