suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# 1. Read the count matrix and the sample sheet
counts <- read.delim("results/counts/count_matrix.tsv", row.names = 1, check.names = FALSE)
meta   <- read.delim("metadata/samples.tsv", stringsAsFactors = FALSE)

# The sample sheet is in a different order than the count columns: reorder it to match
meta <- meta[match(colnames(counts), meta$run_accession), ]
rownames(meta) <- meta$run_accession
stopifnot(identical(rownames(meta), colnames(counts)))
meta$condition <- factor(meta$condition, levels = c("WT", "KO"))   # WT = reference level

# 2. Build the DESeq2 object and drop genes with almost no reads
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes kept:", nrow(dds), "of", nrow(counts), "\n")
cat("Counted pairs per sample:\n"); print(colSums(counts(dds)))

# 3. Variance-stabilising transform (puts counts on a log-like scale for plotting)
vsd <- vst(dds, blind = TRUE)

# 4. PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 4) + geom_text(vjust = -1, size = 3.5) +
  xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw()
ggsave("results/deseq2_qc/PCA.png", p, width = 6, height = 5, dpi = 150)

# 5. Sample-to-sample distances and correlations
d <- dist(t(assay(vsd)))
pheatmap(as.matrix(d), clustering_distance_rows = d, clustering_distance_cols = d,
         filename = "results/deseq2_qc/sample_distance_heatmap.png")
cat("\nSample correlations (VST):\n"); print(round(cor(assay(vsd)), 4))

# PCA numbers on the 500 most variable genes (first 3 PCs)
top <- head(order(apply(assay(vsd), 1, var), decreasing = TRUE), 500)
pc  <- prcomp(t(assay(vsd)[top, ]))
cat("\nVariance explained by PC1-3 (%):\n"); print(round(100 * summary(pc)$importance[2, 1:3], 1))
cat("\nSample coordinates:\n"); print(round(pc$x[, 1:3], 1))
