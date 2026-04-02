#!/usr/bin/env Rscript
# Enrichment plots for EYH and LCZ compounds
# Creates network plots, treeplot, and other visualizations

options(stringsAsFactors = F)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)

if (interactive()) {
  setwd('~/d/sci/src/kd_moa')
}

# Read the ranked data (created by previous GSEA analysis)
eyh = read_tsv('output/tmt_18h_eyh_ranked.tsv', show_col_types = FALSE)
lcz = read_tsv('output/tmt_18h_lcz_ranked.tsv', show_col_types = FALSE)

cat("Data loaded:\n")
cat("EYH:", nrow(eyh), "proteins\n")
cat("LCZ:", nrow(lcz), "proteins\n\n")

# Prepare ranked gene lists (signed by fold-change direction)
# This preserves directionality: negative = downregulated, positive = upregulated
eyh_vec = sign(eyh$l2fc) * eyh$rank_metric
names(eyh_vec) = eyh$gene
eyh_vec = sort(eyh_vec, decreasing = TRUE)

lcz_vec = sign(lcz$l2fc) * lcz$rank_metric
names(lcz_vec) = lcz$gene
lcz_vec = sort(lcz_vec, decreasing = TRUE)

# Remove any NA or empty gene names
eyh_vec = eyh_vec[!is.na(names(eyh_vec)) & names(eyh_vec) != ""]
lcz_vec = lcz_vec[!is.na(names(lcz_vec)) & names(lcz_vec) != ""]

cat("Ranked vectors created:\n")
cat("EYH:", length(eyh_vec), "genes\n")
cat("LCZ:", length(lcz_vec), "genes\n\n")

# Run compareCluster with GSEA on GO Biological Process
cat("Running compareCluster with gseGO...\n")
cat("This may take several minutes...\n\n")

set.seed(42)
edo = compareCluster(
  geneCluster = list(EYH = eyh_vec, LCZ = lcz_vec),
  fun = "gseGO",
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",  # Biological Process
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  minGSSize = 15,
  maxGSSize = 500
)

cat("compareCluster complete!\n")
cat("Significant terms found:\n")
cat("  Total:", nrow(edo@compareClusterResult), "\n")
summary_counts = edo@compareClusterResult %>%
  group_by(Cluster) %>%
  summarize(n = n())
print(summary_counts)

# Save the compareCluster results
write_tsv(as.data.frame(edo), 'output/compareCluster_results.tsv')

# Calculate pairwise term similarity (REQUIRED for emapplot)
cat("\nCalculating pairwise term similarity...\n")
edo = pairwise_termsim(edo)

cat("Similarity calculation complete!\n\n")

# Run individual GSEA for cnetplots
cat("Running individual gseGO for EYH...\n")
set.seed(42)
eyh_gsea_go = gseGO(
  geneList = eyh_vec,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  minGSSize = 15,
  maxGSSize = 500
)
eyh_gsea_go = setReadable(eyh_gsea_go, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")

cat("Running individual gseGO for LCZ...\n")
set.seed(42)
lcz_gsea_go = gseGO(
  geneList = lcz_vec,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "fdr",
  minGSSize = 15,
  maxGSSize = 500
)
lcz_gsea_go = setReadable(lcz_gsea_go, OrgDb = org.Mm.eg.db, keyType = "SYMBOL")

cat("Individual GSEA analyses complete!\n\n")

# ============================================================================
# PLOT 1: Dot plot - Simple overview
# ============================================================================
cat("Creating Plot 1: Dot plot...\n")

resx = 300
png('display_items/enrichplot_dotplot.png', width = 10 * resx, height = 8 * resx, res = resx)

p1 = dotplot(edo, showCategory = 20, font.size = 8) +
  ggtitle("Top 20 Enriched GO Terms per Compound") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p1)
dev.off()

cat("  Saved: display_items/enrichplot_dotplot.png\n")

# ============================================================================
# PLOT 2: cnetplot for EYH
# ============================================================================
cat("Creating Plot 2: cnetplot for EYH...\n")

png('display_items/enrichplot_cnet_eyh.png', width = 14 * resx, height = 14 * resx, res = resx)

p2 = cnetplot(eyh_gsea_go,
              showCategory = 10,
              node_label = "all") +
  ggtitle("EYH - Gene-Concept Network") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right")

print(p2)
dev.off()

cat("  Saved: display_items/enrichplot_cnet_eyh.png\n")

# ============================================================================
# PLOT 3: cnetplot for LCZ
# ============================================================================
cat("Creating Plot 3: cnetplot for LCZ...\n")

png('display_items/enrichplot_cnet_lcz.png', width = 14 * resx, height = 14 * resx, res = resx)

p3 = cnetplot(lcz_gsea_go,
              showCategory = 10,
              node_label = "all") +
  ggtitle("LCZ - Gene-Concept Network") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right")

print(p3)
dev.off()

cat("  Saved: display_items/enrichplot_cnet_lcz.png\n")
