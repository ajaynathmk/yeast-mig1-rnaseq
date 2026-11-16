#!/usr/bin/env Rscript
# DESeq2 differential expression: MIG1 knockout vs wild type
# Usage: Rscript scripts/02_DESeq2.R counts.tsv samples.tsv annotation.gtf outdir padj lfc

# ---- 1. arguments ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript scripts/02_DESeq2.R counts.tsv samples.tsv annotation.gtf outdir padj lfc")
}

counts_file <- args[1]
meta_file   <- args[2]
gtf_file    <- args[3]
out_dir     <- args[4]
padj_cutoff <- as.numeric(args[5])
lfc_cutoff  <- as.numeric(args[6])
min_count   <- 10

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. libraries ----
suppressPackageStartupMessages({
  library(DESeq2)
  library(rtracklayer)
  library(apeglm)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  library(enrichplot)
})
# ---- 3. load data ----

counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
meta   <- read.delim(meta_file)

meta <- meta[match(colnames(counts), meta$run_accession), ]
rownames(meta) <- meta$run_accession

stopifnot(identical(rownames(meta), colnames(counts)))

# ---- 4. build dds and filter ----

meta$condition <- factor(meta$condition, levels = c("WT", "KO"))

dds <- DESeqDataSetFromMatrix(as.matrix(counts), meta, design = ~ condition)

n_min <- min(table(dds$condition))
keep  <- rowSums(counts(dds) >= min_count) >= n_min
dds <- dds[keep, ]


# ---- 5. run DESeq2 ----

dds <- DESeq(dds)

# ---- 6. results and shrinkage ----

res <- results(dds, contrast = c("condition", "KO", "WT"), alpha = padj_cutoff)
res_shrunk <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm", res = res)

# ---- 7. annotate and save tables ----

# gene names and biotypes from the GTF
gtf <- rtracklayer::import(gtf_file)
genes_df <- as.data.frame(gtf[gtf$type == "gene"])

# lookup table: one row per gene, label = name if it exists, otherwise the ID
gene_info <- genes_df[, c("gene_id", "gene_name", "gene_biotype")]
gene_info$label <- ifelse(is.na(gene_info$gene_name), gene_info$gene_id, gene_info$gene_name)

# one result table with both fold changes (shrunken and unshrunken)
res_table <- data.frame(gene_id       = rownames(res),
                        baseMean      = res$baseMean,
                        log2FC_mle    = res$log2FoldChange,
                        log2FC_shrunk = res_shrunk$log2FoldChange,
                        lfcSE         = res_shrunk$lfcSE,
                        pvalue        = res_shrunk$pvalue,
                        padj          = res_shrunk$padj)
res_table <- merge(res_table, gene_info, by = "gene_id", all.x = TRUE)
stopifnot(nrow(res_table) == nrow(res))   # the merge must not change the number of rows

# apply the cutoffs
sig <- !is.na(res_table$padj) & res_table$padj < padj_cutoff &
       abs(res_table$log2FC_shrunk) > lfc_cutoff
res_table$direction <- ifelse(sig & res_table$log2FC_shrunk > 0, "up in KO",
                       ifelse(sig & res_table$log2FC_shrunk < 0, "down in KO", "not significant"))

# sort by significance and save
res_table <- res_table[order(res_table$padj), ]
write.csv(res_table, file.path(out_dir, "deseq2_results_annotated.csv"), row.names = FALSE)
write.csv(as.data.frame(counts(dds, normalized = TRUE)), file.path(out_dir, "normalized_counts.csv"))

cat("up in KO:", sum(res_table$direction == "up in KO"),
    " down in KO:", sum(res_table$direction == "down in KO"), "\n")

# ---- 8. plots ----

# diagnostic plots: dispersion, MA, p-value histogram
png(file.path(out_dir, "diagnostics.png"), width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3))
plotDispEsts(dds)
plotMA(res_shrunk, ylim = c(-5, 5))
hist(res$pvalue[res$baseMean > 1], breaks = 20, col = "grey70",
     main = "p-value histogram", xlab = "p-value")
dev.off()

# volcano plot
volcano <- EnhancedVolcano(res_table,
                           lab = res_table$label,
                           x = "log2FC_shrunk",
                           y = "padj",
                           selectLab = c("MIG1", head(res_table$label, 14)),
                           pCutoff = padj_cutoff,
                           FCcutoff = lfc_cutoff,
                           title = "MIG1 deletion vs wild type")
ggsave(file.path(out_dir, "volcano.png"), volcano, width = 9, height = 8, dpi = 150)

# ---- 9. GO enrichment ----

# genes: all tested genes as background, significant genes split by direction
universe   <- res_table$gene_id[!is.na(res_table$padj)]
up_genes   <- res_table$gene_id[res_table$direction == "up in KO"]
down_genes <- res_table$gene_id[res_table$direction == "down in KO"]

# run enrichment for one ontology, for up and down genes; save the tables
run_go <- function(ont) {
  up   <- enrichGO(gene = up_genes,   universe = universe, OrgDb = org.Sc.sgd.db,
                   keyType = "ORF", ont = ont, pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  down <- enrichGO(gene = down_genes, universe = universe, OrgDb = org.Sc.sgd.db,
                   keyType = "ORF", ont = ont, pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  write.csv(as.data.frame(up),   file.path(out_dir, paste0("GO_", ont, "_up.csv")))
  write.csv(as.data.frame(down), file.path(out_dir, paste0("GO_", ont, "_down.csv")))
  list(up = up, down = down)
}

# dotplot of the top terms, saved to a file (skipped if there are none)
plot_go <- function(x, ont, direction, height = 6) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  p <- dotplot(x, showCategory = 15, title = paste0("GO ", ont, ": ", direction, " in KO"))
  ggsave(file.path(out_dir, paste0("GO_", ont, "_", direction, "_dotplot.png")),
         p, width = 8, height = height, dpi = 300)
}

# run all three ontologies
for (ont in c("BP", "MF", "CC")) {
  go <- run_go(ont)
  plot_go(go$up,   ont, "up")
  plot_go(go$down, ont, "down")
}
