#!/usr/bin/env Rscript
# Sample QC: PCA, sample distances and correlations
# Usage: Rscript scripts/01_sample_qc.R counts.tsv samples.tsv outdir

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript scripts/01_sample_qc.R counts.tsv samples.tsv outdir")
}

counts_file <- args[1]
meta_file   <- args[2]
out_dir     <- args[3]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

wrap_sub <- function(x, width = 70) paste(strwrap(x, width), collapse = "\n")

# 1. Read the count matrix and the sample sheet
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)
meta   <- read.delim(meta_file, stringsAsFactors = FALSE)

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
write.csv(data.frame(sample = colnames(dds), counted_pairs = colSums(counts(dds))),
          file.path(out_dir, "counted_pairs.csv"), row.names = FALSE)

# 3. Variance-stabilising transform (puts counts on a log-like scale for plotting)
vsd <- vst(dds, blind = TRUE)

# 4. PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = c("WT" = "#0072B2", "KO" = "#D55E00"), name = "Genotype") +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(title = "Principal component analysis of samples",
       subtitle = wrap_sub(paste0("Variance-stabilised counts, 500 most variable of ", nrow(vsd),
                                  " genes (total count >= 10)"), width = 55),
       x = paste0("PC1: ", pv[1], "% variance"),
       y = paste0("PC2: ", pv[2], "% variance")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "PCA.png"), p, width = 6.5, height = 5, dpi = 300)
ggsave(file.path(out_dir, "PCA.pdf"), p, width = 6.5, height = 5)

# 5. Sample-to-sample distances and correlations
d <- dist(t(assay(vsd)))
annot <- data.frame(Genotype = meta$condition, row.names = rownames(meta))
annot_cols <- list(Genotype = c(WT = "#0072B2", KO = "#D55E00"))
hm <- pheatmap(as.matrix(d), clustering_distance_rows = d, clustering_distance_cols = d,
               annotation_col = annot, annotation_row = annot, annotation_colors = annot_cols,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
               main = "Sample-to-sample distances (Euclidean, VST counts)",
               fontsize = 10, border_color = NA, silent = TRUE)
png(file.path(out_dir, "sample_distance_heatmap.png"), width = 6.5, height = 5.5, units = "in", res = 300)
grid::grid.draw(hm$gtable); invisible(dev.off())
pdf(file.path(out_dir, "sample_distance_heatmap.pdf"), width = 6.5, height = 5.5)
grid::grid.draw(hm$gtable); invisible(dev.off())
cat("\nSample correlations (VST):\n"); print(round(cor(assay(vsd)), 4))
write.csv(round(cor(assay(vsd)), 4), file.path(out_dir, "sample_correlations.csv"))


# PCA numbers on the 500 most variable genes (first 3 PCs)
top <- head(order(apply(assay(vsd), 1, var), decreasing = TRUE), 500)
pc  <- prcomp(t(assay(vsd)[top, ]))
cat("\nVariance explained by PC1-3 (%):\n"); print(round(100 * summary(pc)$importance[2, 1:3], 1))
cat("\nSample coordinates:\n"); print(round(pc$x[, 1:3], 1))
write.csv(round(pc$x[, 1:3], 2), file.path(out_dir, "pca_coordinates.csv"))