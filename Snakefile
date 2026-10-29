"""RNA-seq pipeline for a yeast gene-deletion strain vs wild type.

FASTQ (ENA) -> FastQC -> STAR -> featureCounts -> MultiQC -> DESeq2 -> GO enrichment
Everything gene-specific lives in config.yaml and metadata/samples.tsv.
"""
import os
import pandas as pd

configfile: "config.yaml"

# ------------------------------------------------------------------
# Sample sheet
# ------------------------------------------------------------------
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("run_accession", drop=False)

valid = {config["reference_level"], config["contrast_level"]}
bad = samples[~samples["condition"].isin(valid)]
if len(bad):
    raise ValueError(
        config["samples"] + ": these runs have a condition outside " + str(sorted(valid))
        + " (CHECK / blank?): " + str(bad.index.tolist())
    )

SAMPLES = samples.index.tolist()

# fastq file name -> download URL (ENA lists URLs without a scheme)
FASTQ_URLS = {}
for _, row in samples.iterrows():
    for u in row["fastq_ftp"].split(";"):
        FASTQ_URLS[os.path.basename(u)] = "https://" + u
FASTQC_BASES = [n[: -len(".fastq.gz")] for n in FASTQ_URLS]


def run_fastqs(wildcards):
    """All FASTQ files of one run, sorted so _1 comes before _2."""
    names = sorted(os.path.basename(u) for u in samples.loc[wildcards.run, "fastq_ftp"].split(";"))
    return ["data/raw/" + n for n in names]


GENOME = "data/reference/genome.fa"
GTF = "data/reference/annotation.gtf"
PAIRED_FLAG = "-p --countReadPairs" if config["paired_end"] else ""


# ------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------
rule all:
    input:
        "results/multiqc/multiqc_report.html",
        "results/deseq2/deseq2_results.csv",
        "results/deseq2/plots/PCA.png",
        "results/enrichment/GO_BP_up.csv",
        "results/enrichment/GO_BP_down.csv",
        "results/enrichment/GSEA_GO_BP.csv",


# ------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------
rule download_fastq:
    output:
        "data/raw/{name}.fastq.gz",
    params:
        url=lambda w: FASTQ_URLS[w.name + ".fastq.gz"],
    log:
        "logs/download/{name}.log",
    shell:
        "wget -q --tries=3 {params.url} -O {output} 2> {log}"


rule download_genome:
    output:
        GENOME,
    params:
        url=config["genome_url"],
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz"


rule download_gtf:
    output:
        GTF,
    params:
        url=config["gtf_url"],
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz"


# ------------------------------------------------------------------
# 2. Quality control of raw reads
# ------------------------------------------------------------------
rule fastqc:
    input:
        "data/raw/{base}.fastq.gz",
    output:
        html="results/qc/fastqc/{base}_fastqc.html",
        zip="results/qc/fastqc/{base}_fastqc.zip",
    threads: 2
    log:
        "logs/fastqc/{base}.log",
    shell:
        "fastqc -t {threads} -o results/qc/fastqc {input} > {log} 2>&1"


# ------------------------------------------------------------------
# 3. Alignment
# ------------------------------------------------------------------
rule star_index:
    input:
        genome=GENOME,
        gtf=GTF,
    output:
        directory("data/reference/star_index"),
    threads: 4
    log:
        "logs/star_index.log",
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate --runThreadN {threads} \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {config[sjdb_overhang]} \
             --genomeSAindexNbases {config[star_sa_index_nbases]} \
             --outFileNamePrefix data/reference/star_index_ > {log} 2>&1
        """


rule star_align:
    input:
        fastq=run_fastqs,
        index="data/reference/star_index",
    output:
        bam="results/star/{run}/Aligned.sortedByCoord.out.bam",
        final_log="results/star/{run}/Log.final.out",
    threads: 4
    log:
        "logs/star/{run}.log",
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesIn {input.fastq} \
             --readFilesCommand zcat \
             --outFileNamePrefix results/star/{wildcards.run}/ \
             --outSAMtype BAM SortedByCoordinate > {log} 2>&1
        """


rule bam_index:
    input:
        "results/star/{run}/Aligned.sortedByCoord.out.bam",
    output:
        "results/star/{run}/Aligned.sortedByCoord.out.bam.bai",
    shell:
        "samtools index {input}"


# ------------------------------------------------------------------
# 4. Counting
# ------------------------------------------------------------------
rule featurecounts:
    input:
        bams=expand("results/star/{run}/Aligned.sortedByCoord.out.bam", run=SAMPLES),
        bais=expand("results/star/{run}/Aligned.sortedByCoord.out.bam.bai", run=SAMPLES),
        gtf=GTF,
    output:
        counts="results/counts/counts.txt",
        summary="results/counts/counts.txt.summary",
    params:
        strand=config["strandedness"],
        paired=PAIRED_FLAG,
    threads: 4
    log:
        "logs/featurecounts.log",
    shell:
        """
        featureCounts -T {threads} -s {params.strand} {params.paired} \
            -t exon -g gene_id -a {input.gtf} \
            -o {output.counts} {input.bams} > {log} 2>&1
        """


# ------------------------------------------------------------------
# 5. Combined QC report
# ------------------------------------------------------------------
rule multiqc:
    input:
        expand("results/qc/fastqc/{base}_fastqc.zip", base=FASTQC_BASES),
        expand("results/star/{run}/Log.final.out", run=SAMPLES),
        "results/counts/counts.txt.summary",
    output:
        "results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        "multiqc results/qc results/star results/counts -o results/multiqc -f > {log} 2>&1"


# ------------------------------------------------------------------
# 6. Statistics
# ------------------------------------------------------------------
rule deseq2:
    input:
        counts="results/counts/counts.txt",
        samples=config["samples"],
    output:
        results="results/deseq2/deseq2_results.csv",
        sig="results/deseq2/deseq2_significant.csv",
        norm="results/deseq2/normalized_counts.csv",
        pca="results/deseq2/plots/PCA.png",
    log:
        "logs/deseq2.log",
    shell:
        """
        Rscript --vanilla scripts/deseq2_analysis.R \
            {input.counts} {input.samples} results/deseq2 \
            {config[target_gene_id]} {config[reference_level]} {config[contrast_level]} \
            {config[padj_cutoff]} {config[lfc_cutoff]} {config[min_count]} > {log} 2>&1
        """


rule enrichment:
    input:
        "results/deseq2/deseq2_results.csv",
    output:
        up="results/enrichment/GO_BP_up.csv",
        down="results/enrichment/GO_BP_down.csv",
        gsea="results/enrichment/GSEA_GO_BP.csv",
    log:
        "logs/enrichment.log",
    shell:
        """
        Rscript --vanilla scripts/enrichment.R \
            {input} results/enrichment {config[padj_cutoff]} {config[lfc_cutoff]} > {log} 2>&1
        """
