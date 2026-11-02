suppressPackageStartupMessages({ library(DESeq2); library(ggplot2) })

# 1. Same loading and ordering as the QC script
counts <- read.delim("results/counts/count_matrix.tsv", row.names = 1, check.names = FALSE)
meta   <- read.delim("metadata/samples.tsv", stringsAsFactors = FALSE)
meta   <- meta[match(colnames(counts), meta$run_accession), ]
rownames(meta) <- meta$run_accession
stopifnot(identical(rownames(meta), colnames(counts)))
meta$condition <- factor(meta$condition, levels = c("WT", "KO"))

dds <- DESeqDataSetFromMatrix(as.matrix(counts), meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 2. Fit the model: normalisation (size factors), per-gene dispersion, Wald test
dds <- DESeq(dds)
cat("Size factors:\n"); print(round(sizeFactors(dds), 3))
cat("Coefficients:", paste(resultsNames(dds), collapse = ", "), "\n")

# 3. Results KO vs WT (positive log2FC = higher in KO), plus shrunken fold changes
res    <- results(dds, contrast = c("condition", "KO", "WT"), alpha = 0.05)
resLFC <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm", res = res)
cat("\nSummary at padj < 0.05 (unshrunken):\n"); summary(res)

out <- as.data.frame(resLFC)
out$gene_id <- rownames(out)
out$log2FC_unshrunk <- res$log2FoldChange
out <- out[order(out$padj), c("gene_id", "baseMean", "log2FoldChange", "log2FC_unshrunk", "lfcSE", "pvalue", "padj")]
write.csv(out, "results/deseq2/deseq2_results.csv", row.names = FALSE)
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "results/deseq2/normalized_counts.csv")

# 4. Headline numbers and the positive control
sig <- subset(out, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nSignificant (padj<0.05, |shrunken log2FC|>1): up in KO =", sum(sig$log2FoldChange > 0),
    ", down in KO =", sum(sig$log2FoldChange < 0), "\n")
cat("\nMIG1 (YGL035C):\n"); print(out[out$gene_id == "YGL035C", ])

# 5. Volcano plot
v <- out[!is.na(out$padj), ]
v$status <- ifelse(v$padj < 0.05 & abs(v$log2FoldChange) > 1, "significant", "not significant")
p <- ggplot(v, aes(log2FoldChange, -log10(padj), colour = status)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c("significant" = "firebrick", "not significant" = "grey60")) +
  labs(x = "log2 fold change (KO vs WT, shrunken)", y = "-log10 adjusted p") + theme_bw()
ggsave("results/deseq2/volcano.png", p, width = 6, height = 5, dpi = 150)
