#!/usr/bin/env Rscript
# Differential expression with DESeq2: deletion strain (KO) vs wild type (WT).
# Args: counts.txt samples.tsv outdir target_gene_id ref_level contrast_level padj lfc min_count

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[1]; samples_file <- args[2]; outdir <- args[3]
target_id <- args[4]; ref_level <- args[5]; contrast_level <- args[6]
padj_cut <- as.numeric(args[7]); lfc_cut <- as.numeric(args[8]); min_count <- as.numeric(args[9])

plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Read featureCounts matrix ------------------------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, 7:ncol(fc)])            # columns 1-6 are Geneid, Chr, Start, End, Strand, Length
rownames(counts) <- fc$Geneid
colnames(counts) <- sub(".*/([EDS]RR[0-9]+)/.*", "\\1", colnames(counts))   # BAM path -> run accession

## ---- 2. Read sample sheet, align order ---------------------------------------
meta <- read.delim(samples_file, stringsAsFactors = FALSE)
rownames(meta) <- meta$run_accession
stopifnot("Count columns not found in sample sheet" = all(colnames(counts) %in% rownames(meta)))
meta <- meta[colnames(counts), ]
meta$condition <- factor(meta$condition, levels = c(ref_level, contrast_level))  # first level = reference
cat("Samples per condition:\n"); print(table(meta$condition))

## ---- 3. DESeq2 -------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= min_count, ]
cat("Genes kept after low-count filter:", nrow(dds), "\n")
dds <- DESeq(dds)                                 # size factors -> dispersions -> Wald test

coef_name <- paste0("condition_", contrast_level, "_vs_", ref_level)
stopifnot(coef_name %in% resultsNames(dds))
res_raw <- results(dds, name = coef_name, alpha = padj_cut)            # unshrunk (for stat / p-values)
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")           # shrunken LFC (for plots / ranking)

## ---- 4. Assemble + annotate result table --------------------------------------
out <- as.data.frame(res_raw)
out$log2FC_shrunk <- res_shr$log2FoldChange[match(rownames(out), rownames(res_shr))]
out$gene_id <- rownames(out)

gene_name <- description <- rep(NA_character_, nrow(out))
try({
  suppressPackageStartupMessages(library(org.Sc.sgd.db))
  gene_name   <- suppressMessages(mapIds(org.Sc.sgd.db, keys = out$gene_id, keytype = "ORF",
                                         column = "GENENAME", multiVals = "first"))
  description <- suppressMessages(mapIds(org.Sc.sgd.db, keys = out$gene_id, keytype = "ORF",
                                         column = "DESCRIPTION", multiVals = "first"))
}, silent = TRUE)
out$gene_name <- unname(gene_name)
out$description <- unname(description)
out$label <- ifelse(is.na(out$gene_name), out$gene_id, out$gene_name)

out$significant <- !is.na(out$padj) & out$padj < padj_cut & abs(out$log2FC_shrunk) >= lfc_cut
out$direction <- ifelse(!out$significant, "NS", ifelse(out$log2FC_shrunk > 0, "Up", "Down"))
out <- out[order(out$padj), c("gene_id", "gene_name", "description", "baseMean", "log2FoldChange",
                              "log2FC_shrunk", "lfcSE", "stat", "pvalue", "padj", "direction", "significant")]

write.csv(out, file.path(outdir, "deseq2_results.csv"), row.names = FALSE)
write.csv(subset(out, significant), file.path(outdir, "deseq2_significant.csv"), row.names = FALSE)

nc <- counts(dds, normalized = TRUE)
write.csv(data.frame(gene_id = rownames(nc), nc, check.names = FALSE),
          file.path(outdir, "normalized_counts.csv"), row.names = FALSE)

cat(sprintf("\nSignificant genes (padj < %g and |shrunk LFC| >= %g): %d up, %d down\n",
            padj_cut, lfc_cut, sum(out$direction == "Up"), sum(out$direction == "Down")))

## ---- 5. Sanity check: the deleted gene should be (almost) gone ---------------
chk <- file.path(outdir, "target_gene_check.txt")
if (target_id %in% rownames(nc)) {
  tg <- out[out$gene_id == target_id, ]
  sink(chk)
  cat("Target gene:", target_id, "\nNormalized counts per sample:\n"); print(round(nc[target_id, ], 1))
  cat("\nlog2FC (unshrunk):", tg$log2FoldChange, " padj:", tg$padj, "\n")
  sink()
  cat(readLines(chk), sep = "\n")
} else {
  cat("WARNING: target gene", target_id, "not in the filtered count matrix (zero reads in every sample?).\n")
}

## ---- 6. QC and exploratory plots ----------------------------------------------
vsd <- vst(dds, blind = TRUE)
cols <- c("#4C78A8", "#E45756"); names(cols) <- c(ref_level, contrast_level)

# PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct <- round(100 * attr(pca, "percentVar"), 1)
p <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 4) + geom_text_repel(show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"), title = "PCA of VST counts") +
  theme_bw()
ggsave(file.path(plot_dir, "PCA.png"), p, width = 6, height = 5, dpi = 300)

# Sample distance heatmap
d <- as.matrix(dist(t(assay(vsd))))
ann <- data.frame(condition = meta$condition, row.names = rownames(meta))
pheatmap(d, annotation_col = ann, annotation_colors = list(condition = cols),
         main = "Sample-to-sample distance", filename = file.path(plot_dir, "sample_distance_heatmap.png"),
         width = 6, height = 5)

# Dispersion estimates
png(file.path(plot_dir, "dispersion.png"), width = 1500, height = 1200, res = 250); plotDispEsts(dds); dev.off()

# p-value histogram (should be flat + a spike near 0 if there is signal)
p <- ggplot(subset(out, baseMean > 1), aes(pvalue)) + geom_histogram(bins = 40, fill = "grey60", colour = "white") +
  labs(title = "p-value distribution", x = "p-value", y = "genes") + theme_bw()
ggsave(file.path(plot_dir, "pvalue_hist.png"), p, width = 6, height = 4, dpi = 300)

## ---- 7. Result plots ------------------------------------------------------------
dcol <- c(Up = "#E45756", Down = "#4C78A8", NS = "grey75")
lab_df <- rbind(head(subset(out, significant & !is.na(padj)), 15), out[out$gene_id == target_id, ])
lab_df$label <- ifelse(is.na(lab_df$gene_name), lab_df$gene_id, lab_df$gene_name)

# MA plot
p <- ggplot(out, aes(baseMean, log2FC_shrunk, colour = direction)) +
  geom_point(alpha = 0.6, size = 1) + scale_x_log10() + scale_colour_manual(values = dcol) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Mean of normalized counts", y = "Shrunken log2 fold change", title = paste(contrast_level, "vs", ref_level)) +
  theme_bw()
ggsave(file.path(plot_dir, "MA_plot.png"), p, width = 7, height = 5, dpi = 300)

# Volcano plot
vd <- out[!is.na(out$padj), ]
p <- ggplot(vd, aes(log2FC_shrunk, -log10(padj), colour = direction)) +
  geom_point(alpha = 0.6, size = 1) + scale_colour_manual(values = dcol) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2) +
  geom_hline(yintercept = -log10(padj_cut), linetype = 2) +
  geom_text_repel(data = lab_df[!is.na(lab_df$padj), ], aes(label = label), colour = "black",
                  size = 3, max.overlaps = 30) +
  labs(x = "Shrunken log2 fold change", y = "-log10 adjusted p-value", title = paste(contrast_level, "vs", ref_level)) +
  theme_bw()
ggsave(file.path(plot_dir, "volcano.png"), p, width = 7, height = 6, dpi = 300)

# Heatmap of top 30 genes (by padj)
top <- head(subset(out, !is.na(padj))$gene_id, 30)
if (length(top) >= 2) {
  mat <- assay(vsd)[top, , drop = FALSE]
  rownames(mat) <- ifelse(is.na(out$gene_name[match(top, out$gene_id)]), top, out$gene_name[match(top, out$gene_id)])
  pheatmap(mat, scale = "row", annotation_col = ann, annotation_colors = list(condition = cols),
           main = "Top 30 genes by adjusted p-value (row z-score of VST)",
           filename = file.path(plot_dir, "top30_heatmap.png"), width = 7, height = 8)
}

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("DESeq2 analysis finished.\n")
