import os
import pandas as pd

configfile: "config.yaml"

samples = pd.read_csv(config["samples"], sep="\t")
SAMPLES = samples["run_accession"].tolist()
READS   = ["1", "2"]

# file name -> download URL, taken from the fastq_ftp column of the sample sheet
FASTQ_URLS = {}
for ftp_list in samples["fastq_ftp"]:
    for u in ftp_list.split(";"):
        FASTQ_URLS[os.path.basename(u)] = u


rule all:
    input:
        "results/multiqc_raw/multiqc_report.html",
        "results/multiqc_trimmed/multiqc_report.html",
        "results/counts/count_matrix.tsv",
        "results/sample_qc/PCA.png",
        "results/deseq2/deseq2_results_annotated.csv"


# ---- 1. downloads ----------------------------------------------------
rule download_fastq:
    output:
        "data/raw/{name}.fastq.gz"
    params:
        url = lambda wc: FASTQ_URLS[wc.name + ".fastq.gz"]
    shell:
        "scripts/download_one.sh {params.url}"

rule download_genome:
    output:
        "data/reference/genome.fa"
    params:
        url = config["genome_url"]
    shell:
        "wget -q -O {output}.gz {params.url} && gunzip -f {output}.gz"

rule download_gtf:
    output:
        "data/reference/annotation.gtf"
    params:
        url = config["gtf_url"]
    shell:
        "wget -q -O {output}.gz {params.url} && gunzip -f {output}.gz"


# ---- 2. quality control and trimming --------------------------------
rule fastqc_raw:
    input:
        "data/raw/{name}.fastq.gz"
    output:
        html = "results/qc/{name}_fastqc.html",
        zip  = "results/qc/{name}_fastqc.zip"
    threads: 2
    shell:
        "fastqc -t {threads} --outdir results/qc {input}"

rule fastp:
    input:
        r1 = "data/raw/{run}_1.fastq.gz",
        r2 = "data/raw/{run}_2.fastq.gz"
    output:
        r1   = "data/trimmed/{run}_1.fastq.gz",
        r2   = "data/trimmed/{run}_2.fastq.gz",
        html = "results/fastp/{run}.html",
        json = "results/fastp/{run}.json"
    log:
        "results/fastp/{run}.log"
    threads: 4
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --detect_adapter_for_pe --length_required 36 \
              --thread {threads} \
              --html {output.html} --json {output.json} \
              > {log} 2>&1
        """

rule fastqc_trimmed:
    input:
        "data/trimmed/{name}.fastq.gz"
    output:
        html = "results/qc_trimmed/{name}_fastqc.html",
        zip  = "results/qc_trimmed/{name}_fastqc.zip"
    threads: 2
    shell:
        "fastqc -t {threads} --outdir results/qc_trimmed {input}"

rule multiqc_raw:
    input:
        expand("results/qc/{run}_{r}_fastqc.zip", run=SAMPLES, r=READS)
    output:
        "results/multiqc_raw/multiqc_report.html"
    shell:
        "multiqc results/qc -o results/multiqc_raw -f"

rule multiqc_trimmed:
    input:
        expand("results/qc_trimmed/{run}_{r}_fastqc.zip", run=SAMPLES, r=READS),
        expand("results/fastp/{run}.json", run=SAMPLES)
    output:
        "results/multiqc_trimmed/multiqc_report.html"
    shell:
        "multiqc results/qc_trimmed results/fastp -o results/multiqc_trimmed -f"


# ---- 3. alignment ----------------------------------------------------
rule star_index:
    input:
        genome = "data/reference/genome.fa",
        gtf    = "data/reference/annotation.gtf"
    output:
        "data/reference/star_index/SAindex"
    params:
        overhang = config["sjdb_overhang"],
        nbases   = config["star_sa_index_nbases"]
    threads: 4
    log:
        "logs/star_index.log"
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} \
             --genomeDir data/reference/star_index \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.overhang} \
             --genomeSAindexNbases {params.nbases} \
             --outFileNamePrefix data/reference/star_index_ > {log} 2>&1
        """

rule star_align:
    input:
        r1    = "data/trimmed/{run}_1.fastq.gz",
        r2    = "data/trimmed/{run}_2.fastq.gz",
        index = "data/reference/star_index/SAindex"
    output:
        bam = "results/star/{run}/Aligned.sortedByCoord.out.bam",
        log = "results/star/{run}/Log.final.out"
    threads: 4
    log:
        "logs/star_{run}.log"
    shell:
        """
        rm -rf results/star/{wildcards.run}/_STARtmp
        STAR --runThreadN {threads} \
             --genomeDir data/reference/star_index \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outFileNamePrefix results/star/{wildcards.run}/ \
             --outSAMtype BAM SortedByCoordinate \
             --alignIntronMax 2000 --alignMatesGapMax 2000 \
             --limitBAMsortRAM 3000000000 > {log} 2>&1
        """


# ---- 4. counting -----------------------------------------------------
rule featurecounts:
    input:
        bams = expand("results/star/{run}/Aligned.sortedByCoord.out.bam", run=SAMPLES),
        gtf  = "data/reference/annotation.gtf"
    output:
        counts  = "results/counts/counts.txt",
        summary = "results/counts/counts.txt.summary"
    params:
        strand = config["strandedness"]
    threads: 4
    log:
        "logs/featurecounts.log"
    shell:
        """
        featureCounts -T {threads} -p --countReadPairs -s {params.strand} \
            -t exon -g gene_id -a {input.gtf} -o {output.counts} {input.bams} > {log} 2>&1
        """

rule count_matrix:
    input:
        "results/counts/counts.txt"
    output:
        "results/counts/count_matrix.tsv"
    shell:
        r"""
        tail -n +2 {input} | cut -f1,7- \
          | sed '1s#results/star/\(SRR[0-9]*\)/Aligned.sortedByCoord.out.bam#\1#g' \
          > {output}
        """


# ---- 5. analysis -----------------------------------------------------
rule sample_qc:
    input:
        counts = "results/counts/count_matrix.tsv",
        meta   = config["samples"]
    output:
        "results/sample_qc/PCA.png",
        "results/sample_qc/sample_distance_heatmap.png",
        "results/sample_qc/counted_pairs.csv",
        "results/sample_qc/sample_correlations.csv",
        "results/sample_qc/pca_coordinates.csv"
    log:
        "logs/sample_qc.log"
    shell:
        "Rscript scripts/01_sample_qc.R {input.counts} {input.meta} results/sample_qc > {log} 2>&1"

rule deseq2:
    input:
        counts = "results/counts/count_matrix.tsv",
        meta   = config["samples"],
        gtf    = config["gtf"]
    output:
        "results/deseq2/deseq2_results_annotated.csv",
        "results/deseq2/normalized_counts.csv",
        "results/deseq2/diagnostics.png",
        "results/deseq2/volcano.png",
        expand("results/deseq2/GO_{ont}_{direction}.csv", ont=["BP", "MF", "CC"], direction=["up", "down"])
    params:
        padj = config["padj_cutoff"],
        lfc  = config["lfc_cutoff"]
    log:
        "logs/deseq2.log"
    shell:
        "Rscript scripts/02_DESeq2.R {input.counts} {input.meta} {input.gtf} results/deseq2 {params.padj} {params.lfc} > {log} 2>&1"