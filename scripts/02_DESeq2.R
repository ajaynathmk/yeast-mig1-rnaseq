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

# shared theme and saver: every plot is written as PNG (report) and PDF (vector)
base_theme <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        legend.position = "bottom", panel.grid.minor = element_blank())

save_plot <- function(p, name, w = 6.5, h = 4.5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
}

# 8a. dispersion estimates
disp <- data.frame(mean = mcols(dds)$baseMean,
                   gene_wise = mcols(dds)$dispGeneEst,
                   fitted = mcols(dds)$dispFit,
                   final = mcols(dds)$dispersion)
disp <- disp[!is.na(disp$gene_wise) & disp$mean > 0, ]
p_disp <- ggplot(disp, aes(mean)) +
  geom_point(aes(y = gene_wise, colour = "Gene-wise estimate"), size = 0.6, alpha = 0.4) +
  geom_point(aes(y = final, colour = "Final (shrunken)"), size = 0.6, alpha = 0.4) +
  geom_line(aes(y = fitted, colour = "Fitted trend"), linewidth = 0.9) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("Gene-wise estimate" = "grey50",
                                 "Final (shrunken)" = "#0072B2",
                                 "Fitted trend" = "#D55E00"), name = NULL) +
  labs(title = "Dispersion estimates",
       subtitle = paste0(nrow(disp), " genes; design ~ condition"),
       x = "Mean of normalised counts", y = "Dispersion") + base_theme
save_plot(p_disp, "dispersion")

# 8b. MA plot (apeglm-shrunken fold changes)
ma <- as.data.frame(res_shrunk)
ma <- ma[!is.na(ma$log2FoldChange), ]
ma$sig <- ifelse(!is.na(ma$padj) & ma$padj < padj_cutoff & abs(ma$log2FoldChange) > lfc_cutoff,
                 "Significant", "Not significant")
p_ma <- ggplot(ma, aes(baseMean, log2FoldChange, colour = sig)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey40") +
  scale_x_log10() +
  scale_colour_manual(values = c("Not significant" = "grey65", "Significant" = "#D55E00"), name = NULL) +
  labs(title = "MA plot, mig1 KO vs wild type",
       subtitle = paste0("apeglm-shrunken log2 fold change; significant = padj < ", padj_cutoff,
                         " and |log2FC| > ", lfc_cutoff),
       x = "Mean of normalised counts", y = "Shrunken log2 fold change (KO / WT)") + base_theme
save_plot(p_ma, "MA_plot")

# 8c. p-value histogram (genes with mean count > 1, as before)
pv <- data.frame(p = res$pvalue[res$baseMean > 1 & !is.na(res$pvalue)])
p_pv <- ggplot(pv, aes(p)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "#0072B2", colour = "white") +
  labs(title = "Distribution of raw p-values",
       subtitle = paste0(nrow(pv), " genes with mean normalised count > 1; a peak near 0 on a flat background is expected"),
       x = "Raw p-value", y = "Number of genes") + base_theme
save_plot(p_pv, "pvalue_histogram")

# 8d. volcano plot
volcano <- EnhancedVolcano(res_table,
                           lab = res_table$label,
                           x = "log2FC_shrunk",
                           y = "padj",
                           selectLab = c("MIG1", head(res_table$label, 14)),
                           pCutoff = padj_cutoff,
                           FCcutoff = lfc_cutoff,
                           title = "MIG1 deletion vs wild type",
                           subtitle = paste0("apeglm-shrunken log2 fold change; cut-offs: padj < ",
                                             padj_cutoff, ", |log2FC| > ", lfc_cutoff),
                           caption = NULL, legendPosition = "bottom")
ggsave(file.path(out_dir, "volcano.png"), volcano, width = 8, height = 7, dpi = 300)
ggsave(file.path(out_dir, "volcano.pdf"), volcano, width = 8, height = 7)


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
go_names <- c(BP = "Biological Process", MF = "Molecular Function", CC = "Cellular Component")

plot_go <- function(x, ont, direction) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(invisible(NULL))
  n_terms <- min(15, nrow(as.data.frame(x)))
  n_genes <- length(if (direction == "up") up_genes else down_genes)
  p <- dotplot(x, showCategory = n_terms, label_format = 40) +
    labs(title = paste0("GO ", go_names[[ont]], ": genes ", direction, " in KO"),
         subtitle = paste0(n_genes, " significant genes; background = ", length(universe),
                           " tested genes; BH-adjusted p < 0.05"),
         x = "Gene ratio (genes in term / significant genes)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          panel.grid.minor = element_blank())
  h <- max(3.5, 0.32 * n_terms + 2)
  ggsave(file.path(out_dir, paste0("GO_", ont, "_", direction, "_dotplot.png")),
         p, width = 7, height = h, dpi = 300)
  ggsave(file.path(out_dir, paste0("GO_", ont, "_", direction, "_dotplot.pdf")),
         p, width = 7, height = h)
}

# run all three ontologies
for (ont in c("BP", "MF", "CC")) {
  go <- run_go(ont)
  plot_go(go$up,   ont, "up")
  plot_go(go$down, ont, "down")
}
