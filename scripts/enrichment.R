#!/usr/bin/env Rscript
# GO enrichment (up / down genes) + GSEA on the ranked Wald statistic.
# Args: deseq2_results.csv outdir padj lfc

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Sc.sgd.db)
  library(enrichplot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
res_file <- args[1]; outdir <- args[2]
padj_cut <- as.numeric(args[3]); lfc_cut <- as.numeric(args[4])
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

res <- read.csv(res_file, stringsAsFactors = FALSE)
universe <- res$gene_id[!is.na(res$padj)]          # background = all genes that were actually tested
up   <- res$gene_id[res$direction == "Up"]
down <- res$gene_id[res$direction == "Down"]
cat("Up:", length(up), " Down:", length(down), " Universe:", length(universe), "\n")

go_enrich <- function(genes, label) {
  csv <- file.path(outdir, paste0("GO_BP_", label, ".csv"))
  if (length(genes) < 5) {
    cat("Too few", label, "genes for enrichment; writing empty table.\n")
    write.csv(data.frame(), csv, row.names = FALSE); return(invisible(NULL))
  }
  ego <- enrichGO(gene = genes, universe = universe, OrgDb = org.Sc.sgd.db, keyType = "ORF",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  df <- as.data.frame(ego)
  write.csv(df, csv, row.names = FALSE)
  if (nrow(df) > 0) {
    ggsave(file.path(outdir, paste0("GO_BP_", label, "_dotplot.png")),
           dotplot(ego, showCategory = 15) + ggtitle(paste("GO BP -", label, "genes")),
           width = 8, height = 7, dpi = 300)
  }
}
go_enrich(up, "up")
go_enrich(down, "down")

# GSEA: uses ALL genes ranked by Wald statistic, so it does not depend on a cutoff
gsea_csv <- file.path(outdir, "GSEA_GO_BP.csv")
tryCatch({
  ranks <- res$stat; names(ranks) <- res$gene_id
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  gse <- gseGO(geneList = ranks, OrgDb = org.Sc.sgd.db, keyType = "ORF", ont = "BP",
               minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  write.csv(as.data.frame(gse), gsea_csv, row.names = FALSE)
  if (nrow(as.data.frame(gse)) > 0) {
    ggsave(file.path(outdir, "GSEA_GO_BP_dotplot.png"),
           dotplot(gse, showCategory = 12, split = ".sign") + facet_grid(. ~ .sign),
           width = 10, height = 7, dpi = 300)
  }
}, error = function(e) {
  cat("GSEA failed:", conditionMessage(e), "\n"); write.csv(data.frame(), gsea_csv, row.names = FALSE)
})

# KEGG (needs internet; skipped quietly if unavailable). Yeast KEGG gene IDs are the systematic names.
for (lab in c("up", "down")) {
  genes <- get(lab)
  if (length(genes) < 5) next
  tryCatch({
    ek <- enrichKEGG(gene = genes, organism = "sce", keyType = "kegg", universe = universe, pvalueCutoff = 0.05)
    write.csv(as.data.frame(ek), file.path(outdir, paste0("KEGG_", lab, ".csv")), row.names = FALSE)
  }, error = function(e) cat("KEGG", lab, "skipped:", conditionMessage(e), "\n"))
}

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("Enrichment finished.\n")
