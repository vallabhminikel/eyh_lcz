#!/usr/bin/env Rscript
# GSEA analysis of TMT proteomics data for EYH and LCZ compounds
# This script performs gene set enrichment analysis and generates summary reports

options(stringsAsFactors = F)
library(tidyverse)
library(fgsea)
library(msigdbr)

if (interactive()) {
  setwd('~/d/sci/src/kd_moa')
}

msigdbr_collections()

# Helper function for FDR correction and ranking
prepare_gsea_data = function(tmt_data, compound_name) {
  tmt_data %>%
    mutate(
      padj = p.adjust(p_ebm, method = "fdr"),
      rank_metric = abs(l2fc) * -log10(padj + 1e-300),
      compound = compound_name
    ) %>%
    arrange(desc(rank_metric)) -> result

  return(result)
}

# Read TMT data
eyh_raw = read_tsv('output/tmt_18h_eyh.tsv', show_col_types = FALSE)
lcz_raw = read_tsv('output/tmt_18h_lcz.tsv', show_col_types = FALSE)

# Prepare data with FDR correction and ranking
eyh = prepare_gsea_data(eyh_raw, "EYH")
lcz = prepare_gsea_data(lcz_raw, "LCZ")

# Combine for comparison
combined = bind_rows(eyh, lcz)

# Save ranked lists
write_tsv(eyh, 'output/tmt_18h_eyh_ranked.tsv')
write_tsv(lcz, 'output/tmt_18h_lcz_ranked.tsv')

cat("FDR correction summary:\n")
cat("EYH: ", sum(eyh$padj < 0.05), " genes at FDR < 0.05\n")
cat("LCZ: ", sum(lcz$padj < 0.05), " genes at FDR < 0.05\n")
cat("EYH: ", sum(eyh$padj < 0.1), " genes at FDR < 0.1\n")
cat("LCZ: ", sum(lcz$padj < 0.1), " genes at FDR < 0.1\n\n")

hallmark_sets <- msigdbr(species = "Mus musculus", collection = "H")
gobp_sets     <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
kegg_sets     <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:KEGG_MEDICUS")

# Combine gene sets
all_sets = bind_rows(
  hallmark_sets,
  gobp_sets,
  kegg_sets
) %>%
  select(gs_name, gene_symbol)

# Convert to list format for fgsea
pathways = split(all_sets$gene_symbol, all_sets$gs_name)

# Function to run GSEA
run_gsea_analysis = function(data, compound_name) {
  # Create ranked gene list based on signed rank metric
  # Use sign of l2fc to preserve direction
  gene_ranks = setNames(
    sign(data$l2fc) * data$rank_metric,
    data$gene
  )

  # Remove any NAs, empty names, or infinite values
  gene_ranks = gene_ranks[!is.na(names(gene_ranks)) &
                          names(gene_ranks) != "" &
                          is.finite(gene_ranks)]

  # Run fgsea (using fgseaMultilevel, the modern recommended method)
  set.seed(42)  # for reproducibility
  fgsea_results = fgsea(
    pathways = pathways,
    stats = gene_ranks,
    minSize = 15,
    maxSize = 500,
    nPermSimple = 10000  # Increase permutations for unbalanced pathways
  )

  fgsea_results %>%
    as_tibble() %>%
    filter(!is.na(padj)) %>%  # Remove pathways with NA p-values
    arrange(padj, pval) %>%
    mutate(compound = compound_name) %>%
    select(compound, pathway, pval, padj, ES, NES, size, leadingEdge) -> result

  return(result)
}

# Run GSEA for both compounds
cat("Running GSEA for EYH...\n")
eyh_gsea = run_gsea_analysis(eyh, "EYH")

cat("Running GSEA for LCZ...\n")
lcz_gsea = run_gsea_analysis(lcz, "LCZ")

# Combine results
all_gsea = bind_rows(eyh_gsea, lcz_gsea)

# Save GSEA results
write_tsv(eyh_gsea, 'output/gsea_eyh.tsv')
write_tsv(lcz_gsea, 'output/gsea_lcz.tsv')
write_tsv(all_gsea, 'output/gsea_combined.tsv')

# Save top 10 pathways for each compound and combined
eyh_gsea %>%
  head(10) %>%
  write_tsv('output/gsea_eyh_top10.tsv')

lcz_gsea %>%
  head(10) %>%
  write_tsv('output/gsea_lcz_top10.tsv')

all_gsea %>%
  arrange(padj, pval) %>%
  head(10) %>%
  write_tsv('output/gsea_combined_top10.tsv')

cat("\nGSEA summary:\n")
cat("EYH: ", sum(eyh_gsea$padj < 0.05, na.rm=TRUE), " pathways at FDR < 0.05\n")
cat("LCZ: ", sum(lcz_gsea$padj < 0.05, na.rm=TRUE), " pathways at FDR < 0.05\n")
cat("EYH: ", sum(eyh_gsea$padj < 0.25, na.rm=TRUE), " pathways at FDR < 0.25\n")
cat("LCZ: ", sum(lcz_gsea$padj < 0.25, na.rm=TRUE), " pathways at FDR < 0.25\n\n")

# Generate concise report
sink('output/gsea_report.txt')

cat("========================================\n")
cat("GSEA Analysis Report\n")
cat("========================================\n\n")

cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

cat("Analysis Overview:\n")
cat("------------------\n")
cat("Compounds analyzed: EYH, LCZ\n")
cat("Input data: TMT proteomics at 18h timepoint\n")
cat("Gene sets: MSigDB Hallmark, GO:BP, KEGG pathways\n")
cat("Ranking metric: |log2FC| * -log10(padj)\n\n")

cat("Differentially Expressed Proteins (FDR < 0.05):\n")
cat("------------------------------------------------\n")
cat(sprintf("EYH: %d proteins (%.1f%% of %d total)\n",
    sum(eyh$padj < 0.05),
    100*sum(eyh$padj < 0.05)/nrow(eyh),
    nrow(eyh)))
cat(sprintf("LCZ: %d proteins (%.1f%% of %d total)\n\n",
    sum(lcz$padj < 0.05),
    100*sum(lcz$padj < 0.05)/nrow(lcz),
    nrow(lcz)))

cat("Top 10 Proteins by Ranking Metric:\n")
cat("-----------------------------------\n")
cat("\nEYH:\n")
eyh %>%
  select(gene, l2fc, p_ebm, padj, rank_metric) %>%
  head(10) %>%
  mutate(across(where(is.numeric), ~formatC(., format="g", digits=3))) %>%
  capture.output() %>%
  cat(sep="\n")

cat("\n\nLCZ:\n")
lcz %>%
  select(gene, l2fc, p_ebm, padj, rank_metric) %>%
  head(10) %>%
  mutate(across(where(is.numeric), ~formatC(., format="g", digits=3))) %>%
  capture.output() %>%
  cat(sep="\n")

cat("\n\nEnriched Pathways (FDR < 0.05):\n")
cat("--------------------------------\n")

cat("\nEYH - Top 10 pathways:\n")
if (sum(eyh_gsea$padj < 0.05) > 0) {
  eyh_gsea %>%
    filter(padj < 0.05) %>%
    select(pathway, NES, padj, size) %>%
    head(10) %>%
    mutate(across(where(is.numeric), ~formatC(., format="g", digits=3))) %>%
    capture.output() %>%
    cat(sep="\n")
} else {
  cat("  No pathways reached FDR < 0.05\n")
}

cat("\n\nLCZ - Top 10 pathways:\n")
if (sum(lcz_gsea$padj < 0.05) > 0) {
  lcz_gsea %>%
    filter(padj < 0.05) %>%
    select(pathway, NES, padj, size) %>%
    mutate(across(where(is.numeric), ~formatC(., format="g", digits=3))) %>%
    capture.output() %>%
    cat(sep="\n")
} else {
  cat("  No pathways reached FDR < 0.05\n")
}

cat("\n\nPathways Enriched in Both Compounds (FDR < 0.05):\n")
cat("--------------------------------------------------\n")
shared_pathways = inner_join(
  eyh_gsea %>% filter(padj < 0.05) %>% select(pathway, NES_eyh=NES, padj_eyh=padj),
  lcz_gsea %>% filter(padj < 0.05) %>% select(pathway, NES_lcz=NES, padj_lcz=padj),
  by = "pathway"
) %>%
  mutate(concordant = sign(NES_eyh) == sign(NES_lcz)) %>%
  arrange(padj_eyh + padj_lcz)

if (nrow(shared_pathways) > 0) {
  shared_pathways %>%
    head(15) %>%
    mutate(across(where(is.numeric), ~formatC(., format="g", digits=3))) %>%
    capture.output() %>%
    cat(sep="\n")
} else {
  cat("  No shared pathways at FDR < 0.25\n")
}

cat("\n")
sink()

# Generate methods statement
sink('output/gsea_methods.txt')

cat("METHODS STATEMENT\n")
cat("=================\n\n")

cat("Gene Set Enrichment Analysis:\n")
cat("Tandem mass tag (TMT) proteomics data from 18-hour compound treatments were analyzed ")
cat("for differential protein expression. Raw p-values from Empirical Brown's Method (combining ")
cat("multiple peptide-level tests per gene) were adjusted for multiple testing using ")
cat("False discovery rate (FDR) correction. Proteins were ranked by the metric ")
cat("|log2(fold-change)| × -log10(FDR-adjusted p-value), which prioritizes proteins with both ")
cat("large effect sizes and high statistical significance. Gene set enrichment analysis (GSEA) ")
cat("was performed using the fgsea R package (Korotkevich et al., 2021) with the fgseaMultilevel ")
cat("algorithm for adaptive multilevel Monte Carlo sampling. ")
cat("Gene sets were obtained from the Molecular Signatures Database (MSigDB) v2024.1, including ")
cat("Hallmark pathways, Gene Ontology Biological Processes (GO:BP), and KEGG pathways. ")
cat("Pathways with 15-500 genes were tested. GSEA uses a signed ranking metric (preserving ")
cat("fold-change direction) to identify pathways enriched in upregulated or downregulated proteins. ")
cat("Pathway enrichment was considered significant at FDR < 0.05.\n\n")

cat("References:\n")
cat("Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A (2021). ")
cat("Fast gene set enrichment analysis. bioRxiv. doi:10.1101/060012\n")

cat("\n")
sink()

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")
cat("\nOutput files created:\n")
cat("  - output/tmt_18h_eyh_ranked.tsv (EYH protein rankings with padj)\n")
cat("  - output/tmt_18h_lcz_ranked.tsv (LCZ protein rankings with padj)\n")
cat("  - output/gsea_eyh.tsv (EYH GSEA results - all pathways)\n")
cat("  - output/gsea_lcz.tsv (LCZ GSEA results - all pathways)\n")
cat("  - output/gsea_combined.tsv (Combined GSEA results - all pathways)\n")
cat("  - output/gsea_eyh_top10.tsv (Top 10 pathways for EYH)\n")
cat("  - output/gsea_lcz_top10.tsv (Top 10 pathways for LCZ)\n")
cat("  - output/gsea_combined_top10.tsv (Top 10 pathways overall)\n")
cat("  - output/gsea_report.txt (Concise summary report)\n")
cat("  - output/gsea_methods.txt (Methods statement)\n\n")

cat("View the report:\n")
cat("  cat output/gsea_report.txt\n\n")
