# Yeast MIG1 deletion RNA-seq: a reproducible Snakemake pipeline

Differential expression analysis of a *mig1Δ* knockout against wild-type *Saccharomyces cerevisiae*, from raw
sequencing reads to GO enrichment. The whole analysis runs as a Snakemake workflow, and the exact software
environment is available as a Docker image.

This project is a learning project: every step was first run by hand, checked, and only then turned into a
script and a workflow rule. The records of the manual work (strandedness test, investigation of a problem
sample) are kept in [`exploration/`](exploration/).

## Contents

1. [Dataset](#dataset)
2. [Workflow overview](#workflow-overview)
3. [How to run it](#how-to-run-it)
4. [Methods and parameters](#methods-and-parameters)
5. [Quality control findings](#quality-control-findings)
6. [Results](#results)
7. [Limitations](#limitations)
8. [Reproducibility checks](#reproducibility-checks)
9. [Repository layout](#repository-layout)
10. [Troubleshooting](#troubleshooting)
11. [Planned work](#planned-work)
12. [Tools and references](#tools-and-references)

---

## Dataset

| | |
|---|---|
| Source | GEO series **GSE104356** (BioProject **PRJNA412405**) |
| Organism | *Saccharomyces cerevisiae* |
| Samples used | 2 wild type (strain ScY01a): `SRR6108775`, `SRR6108776` |
| | 2 *mig1* knockout: `SRR6108783`, `SRR6108784` |
| Library | paired-end, 151 bp, unstranded |
| Reference | Ensembl Genomes fungi, release 56, assembly **R64-1-1** (FASTA and GTF) |

MIG1 (systematic name `YGL035C`) is the main glucose-repression transcription factor in yeast, so the
experiment has a built-in positive control: MIG1 itself must disappear in the knockout, and known
glucose-repressed genes should change.

The sample sheet is [`metadata/samples.tsv`](metadata/samples.tsv). It was built from the ENA file report for
the project and contains the download URLs of all FASTQ files.

---

## Workflow overview

```
ENA download ──► FastQC ──► fastp trimming ──► FastQC + MultiQC (before/after)
                                  │
                                  ▼
          STAR index (genome + GTF) ──► STAR alignment (per sample)
                                  │
                                  ▼
                featureCounts (gene-level, read pairs)
                                  │
                                  ▼
                       count matrix (genes × samples)
                         │                    │
                         ▼                    ▼
                  sample QC (PCA,        DESeq2 differential
                  distances)             expression + apeglm shrinkage
                                                │
                                                ▼
                                    GO enrichment (BP, MF, CC;
                                    up- and down-regulated genes separately)
```

Each box is a Snakemake rule in [`Snakefile`](Snakefile). The rules are linked only by file names, so
Snakemake works out the order itself and re-runs only what is missing or out of date.

| Rule | What it does | Tool |
|---|---|---|
| `download_fastq` | download one FASTQ file, trying https, http, then ftp, and verify its size against ENA | `scripts/download_one.sh` (wget) |
| `download_genome`, `download_gtf` | fetch the reference | wget |
| `fastqc_raw`, `fastqc_trimmed` | per-file quality reports | FastQC |
| `fastp` | adapter removal and trimming | fastp |
| `multiqc_raw`, `multiqc_trimmed` | combined quality reports | MultiQC |
| `star_index` | build the genome index | STAR |
| `star_align` | align each sample, sorted BAM | STAR |
| `featurecounts` | count read pairs per gene | featureCounts |
| `count_matrix` | tidy the counts into a genes × samples table | shell |
| `sample_qc` | PCA, sample distances, correlations | `scripts/01_sample_qc.R` |
| `deseq2` | differential expression, annotation, plots, GO enrichment | `scripts/02_DESeq2.R` |

---

## How to run it

You need about **15 GB of free disk space** (raw reads, trimmed reads, genome index and BAM files) and a
few GB of RAM. The first run downloads about 5 GB of sequencing data.

### Option 1: Docker

The image contains the tools and the workflow code, not the data. Data and results are mounted from your
machine, so they stay on your disk.

```bash
git clone https://github.com/ajaynathmk/yeast-mig1-rnaseq.git
cd yeast-mig1-rnaseq
docker build -t yeast-mig1-rnaseq .

mkdir -p data results
docker run --rm --user $(id -u):$(id -g) -e HOME=/tmp \
  -v $PWD/data:/pipeline/data -v $PWD/results:/pipeline/results \
  yeast-mig1-rnaseq --cores 4
```

- `--user $(id -u):$(id -g)` makes the container write files as you, so the results belong to your user.
- Anything after the image name goes to Snakemake. For a dry run that lists the jobs without running them, use
  `yeast-mig1-rnaseq -n`.

### Option 2: micromamba

```bash
git clone https://github.com/ajaynathmk/yeast-mig1-rnaseq.git
cd yeast-mig1-rnaseq
micromamba create -n yeast-rnaseq -f environment.lock.yml
micromamba activate yeast-rnaseq
snakemake --cores 4
```

Useful Snakemake commands:

```bash
snakemake -n                  # dry run: which jobs would run
snakemake -n --forceall       # dry run as if nothing existed yet (shows the full 42-job plan)
snakemake --cores 4           # run everything that is missing
snakemake --cores 4 -R deseq2 # force one rule (and what depends on it) to run again
```

### Environment files

| File | Purpose |
|---|---|
| `environment.yml` | the packages the project asks for, without versions; short and readable |
| `environment.lock.yml` | every package pinned to an exact version; this is what the Docker image installs |

### Settings

Values that someone may want to change are in [`config.yaml`](config.yaml): sample sheet, reference URLs,
STAR index settings, strandedness, and the significance cut-offs.

---

## Methods and parameters

**Trimming (fastp).** Adapters detected automatically for paired-end reads (`--detect_adapter_for_pe`);
reads shorter than 36 bases after trimming are removed. Duplicates are **not** removed, as is standard for
RNA-seq (it would flatten the dynamic range of highly expressed genes).

**Alignment (STAR 2.7.10b).** Index built with the annotation (`--sjdbOverhang 150`, read length minus one;
`--genomeSAindexNbases 10`, suited to the 12 Mb yeast genome). Alignment with `--alignIntronMax 2000` and
`--alignMatesGapMax 2000`, because yeast introns are short, and `--limitBAMsortRAM 3000000000` for the BAM
sorting step.

**Counting (featureCounts 2.1.1).** `-p --countReadPairs -s 0 -t exon -g gene_id`: read pairs are counted
once per fragment, the library is treated as **unstranded**, exons are summed per gene. Strandedness was
determined by experiment, not assumed: see the QC findings.

**Differential expression (DESeq2).**
- Design `~ condition` with wild type as the reference level.
- Genes are kept if they have at least 10 reads in at least as many samples as the smallest group (here 2).
  This leaves **5,969** genes. A simpler rule (total of 10 reads over all samples) would keep 6,291, including
  genes with reads in a single sample.
- Design `~ replicate + condition` was considered and rejected: with four samples it leaves one residual
  degree of freedom, and the replicate labels are arbitrary.
- Wald test, Benjamini–Hochberg adjustment.
- Fold changes are shrunk with **apeglm**. A gene is called significant if **padj < 0.05 and |shrunken log2
  fold change| > 1**.
- The result table also contains the unshrunken fold change, the gene name and the gene biotype from the GTF.

**GO enrichment (clusterProfiler).** `enrichGO` with `org.Sc.sgd.db`, gene identifiers as `ORF`, Benjamini–
Hochberg adjustment (cut-off 0.05). The background ("universe") is every gene that was tested (non-missing
adjusted p-value), not the whole genome. Up- and down-regulated genes are tested separately, for the three
ontologies Biological Process, Molecular Function and Cellular Component.

---

## Quality control findings

These are the decisions the data forced, and they are why the pipeline looks the way it does.

**1. The library is unstranded.** Counting the same BAM files with `-s 0`, `-s 1` and `-s 2` assigned 81–84%
of read pairs with `-s 0` and about 46% with either stranded setting. The three settings were compared before
any counting was done for the analysis (see [`exploration/strandedness_test/`](exploration/strandedness_test/)).

**2. Sample SRR6108784 (knockout replicate 2) is of lower quality.**
- Only 36.5% of its reads map uniquely, and 61.7% are discarded by STAR as "too short".
- The unmapped reads contain concatemers of the SMART primer sequence (about 20% of unmapped pairs carry the
  primer motif) and reads that BLAST to *Drosophila melanogaster* (a sample of 10 reads), so the library is
  contaminated with foreign RNA. The GC content of the unmapped reads is about 49%, unlike yeast.
- It yields about 3.0 million counted read pairs against 6.3–7.0 million for the other samples, and a
  DESeq2 size factor of 0.53 against 1.23–1.27.
- In the PCA, this sample separates from the other three on one axis, with the other knockout sample
  (`SRR6108783`) on PC1.
- The sample was **kept in the analysis**, since only four samples exist. All 389 significant genes have both
  knockout samples on the same side of both wild-type samples. This is a necessary check, not a proof; under
  no real difference, the chance of that ordering for one gene is one in three.
- The investigation is documented in
  [`exploration/sample_784_diagnosis/`](exploration/sample_784_diagnosis/).

**3. DESeq2 diagnostics.** The dispersion fit, the MA plot and the p-value histogram are written to
`results/deseq2/diagnostics.png`. 28 genes were flagged as dispersion outliers, mostly heat-shock chaperones
and metabolic genes.

---

## Results

![Volcano plot](results/deseq2/volcano.png)

![PCA](results/sample_qc/PCA.png)

**Overview.** Deleting MIG1 changes the expression of **389 of 5,969 tested genes** (193 up, 196 down;
padj < 0.05, |log2FC| > 1). Without the fold-change requirement, 639 genes pass padj < 0.05 (299 up, 340 down).

**Positive controls.**

| Gene | Observation |
|---|---|
| *MIG1* (`YGL035C`) | almost absent in the knockout: log2 fold change about −7.2 (about 150-fold lower), padj 2e-16 |
| *HXT2*, *HXK1* | glucose-repressed genes, up in the knockout (log2FC +2.7 and +4.0) |
| *SUC2* | modestly up (log2FC +0.8, padj 0.008, below the 2-fold cut-off) |
| *GAL* genes | unchanged |

**Strongest changes.** FIT1, INO1, FIT3, MPC3, FET3, SIT1, PHO84, IMD2, IZH4 and CTR3 are among the top
genes. Of the 389 significant genes, 370 are protein-coding, 12 snoRNA, 3 ncRNA, 2 transposable elements,
1 pseudogene and 1 tRNA.

**GO enrichment** (Biological Process, Molecular Function and Cellular Component agree):
- **Up in the knockout:** metal ion (iron, copper, zinc) transport, hexose and sugar transport, carbon
  metabolism, cell wall and vacuole genes, and the pheromone and mating response.
- **Down in the knockout:** the endoplasmic reticulum protein-folding machinery (ER stress, unfolded protein
  binding, disulfide oxidoreductases, peroxidases, ER lumen).

These are patterns in the data. This analysis does not show *why* the knockout has a metal-starvation-like
signature.

**Result files** (in `results/`):

| File | Content |
|---|---|
| `counts/count_matrix.tsv` | raw gene counts, genes × samples |
| `sample_qc/` | PCA, distance heatmap, correlations, counted pairs per sample |
| `deseq2/deseq2_results_annotated.csv` | all tested genes: fold changes, p-values, name, biotype, direction |
| `deseq2/normalized_counts.csv` | size-factor-normalised counts |
| `deseq2/volcano.png`, `diagnostics.png` | plots |
| `deseq2/GO_{BP,MF,CC}_{up,down}.csv` | enrichment tables, plus dot plots where terms were found |
| `multiqc_raw/`, `multiqc_trimmed/` | combined read quality reports |

---

## Limitations

- **Two replicates per group.** The test has little power, and the results should be read as
  hypothesis-generating.
- **Sample SRR6108784** is contaminated and has about half the usable reads of the others (see above). It
  agrees with the other knockout sample on every significant gene, but a leave-one-out analysis has not yet
  been done.
- **The pheromone and mating up-regulation is unexplained.** Both strains appear to be MATa (a check of the
  mating-type locus and of the alpha-specific genes, which are unchanged), so a mating-type difference does
  not explain it.
- **Growth conditions and the genetic background of the two strains** have not been confirmed against the
  publication.
- **snoRNA and tRNA results** (about 13 genes) are unreliable, since the library type captures them poorly.
- Whether the two replicates of each group are independent biological replicates has not been established.

---

## Reproducibility checks

| Check | Outcome |
|---|---|
| Fresh clone of the repository, whole workflow from the download step | count matrix identical (same numbers; column order differs because it follows the sample sheet) |
| Same run, DESeq2 results | same genes, same direction labels, same GO tables; numerical differences up to about 3e-6 (convergence tolerance of the model fit and of apeglm) |
| Docker image, QC and DESeq2 steps from the count matrix | identical to the manual analysis (maximum absolute difference 0) |

---

## Repository layout

```
├── Snakefile                  the workflow
├── config.yaml                settings (sample sheet, URLs, cut-offs)
├── Dockerfile                 image with the pinned environment and the workflow
├── environment.yml            requested packages
├── environment.lock.yml       exact versions
├── metadata/                  sample sheet and ENA file lists
├── scripts/
│   ├── download_one.sh        robust single-file download (https → http → ftp, size check)
│   ├── 01_sample_qc.R         PCA, sample distances, correlations
│   └── 02_DESeq2.R            differential expression, annotation, plots, GO enrichment
├── notebooks/                 the interactive DESeq2 notebook this script was built from
├── exploration/               manual checks that are not part of the pipeline
│   ├── strandedness_test/     featureCounts with -s 0, 1 and 2
│   ├── sample_784_diagnosis/  STAR runs on the problem sample
│   ├── manual_steps/          the pipeline as separate shell scripts (before Snakemake)
│   └── legacy_download/       earlier download scripts
├── docs/
│   ├── Notes.md               running log of the work
│   └── report/                LaTeX report
├── figures/                   images used in this README
└── results/                   pipeline outputs
```

Not in the repository (listed in `.gitignore`): `data/` (raw and trimmed reads, reference, STAR index),
BAM files, and Snakemake's working folders.

---

## Troubleshooting

- **STAR stops during BAM sorting with a message about `--limitBAMsortRAM`.** The sort needs more memory than
  STAR assumed. The workflow passes `--limitBAMsortRAM 3000000000` (3 GB); raise it if your machine allows.
  A leftover `_STARtmp` folder also stops STAR from starting; the rule removes it before each run.
- **A download fails or stalls.** ENA serves some files only over some protocols. `scripts/download_one.sh`
  tries https, http and ftp in turn and resumes partial files; it counts a file as complete only when its
  size equals the size ENA reports (`metadata/fastq_sizes.tsv`). Snakemake deletes the partial output of a
  failed job; to keep partial files for resuming, run with `--keep-incomplete`.
- **`Permission denied` when running in Docker.** The container user is not your user. Use
  `--user $(id -u):$(id -g) -e HOME=/tmp` as shown above.
- **`micromamba create` reports conflicts or `repo.anaconda.com`.** The environment files use only
  `conda-forge` and `bioconda` (`nodefaults`). Do not add Anaconda's default channel to your configuration.
- **Snakemake says "Nothing to be done" although you changed an R script.** Snakemake only watches the
  files that a rule lists as inputs. Force the rule with `-R sample_qc` or `-R deseq2`.

---

## Planned work

- Compare the results with the gene list deposited by the authors with the GEO series, as a script and
  workflow rule, so that the comparison is reproducible.
- Leave-one-out analysis without `SRR6108784`, and a model with a replicate term, as a sensitivity check.
- Gene set enrichment on the ranked gene list, and KEGG pathways.
- Check growth conditions and strain background against the publication.

---

## Tools and references

Snakemake 7.32.4, STAR 2.7.10b, featureCounts (Subread) 2.1.1, fastp, FastQC, MultiQC, DESeq2, apeglm,
clusterProfiler, enrichplot, EnhancedVolcano, and the *S. cerevisiae* annotation package `org.Sc.sgd.db`.
Exact versions of everything are in [`environment.lock.yml`](environment.lock.yml).

- Mölder F, et al. Sustainable data analysis with Snakemake. *F1000Research* 2021.
- Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 2013.
- Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads
  to genomic features. *Bioinformatics* 2014.
- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with
  DESeq2. *Genome Biology* 2014.
- Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise
  and preserving large differences. *Bioinformatics* 2019.
- Chen S, et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 2018.
- Ewels P, et al. MultiQC: summarize analysis results for multiple tools and samples in a single report.
  *Bioinformatics* 2016.
- Wu T, et al. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *The
  Innovation* 2021.
- The publication associated with GEO series GSE104356 (linked from the series page): TO ADD.