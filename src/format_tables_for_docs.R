# Format GSEA top10 tables for Google Docs
library(tidyverse)

# Function to clean pathway names
clean_pathway <- function(pathway) {
  pathway %>%
    str_remove("^GOBP_") %>%
    str_remove("^HALLMARK_") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

# Function to format scientific notation
format_pval <- function(pval) {
  ifelse(pval < 0.001,
         sprintf("%.2e", pval),
         sprintf("%.4f", pval))
}

# Read and format each file
files <- c("gsea_combined_top10", "gsea_eyh_top10", "gsea_lcz_top10")

for (file in files) {
  cat("\n===========================================\n")
  cat(toupper(str_remove(file, "_top10")), "\n")
  cat("===========================================\n\n")

  df <- read_tsv(paste0("output/", file, ".tsv"), show_col_types = FALSE)

  # Format the table
  formatted <- df %>%
    mutate(
      Pathway = clean_pathway(pathway),
      `P-value` = format_pval(pval),
      `Adj. P-value` = format_pval(padj),
      `Enrich. Score` = round(ES, 3),
      `Norm. Enrich. Score` = round(NES, 3),
      `Gene Set Size` = size
    ) %>%
    select(Compound = compound, Pathway, `P-value`, `Adj. P-value`,
           `Enrich. Score`, `Norm. Enrich. Score`, `Gene Set Size`)

  # Print as tab-separated for easy copying
  write.table(formatted, sep="\t", quote=FALSE, row.names=FALSE)

  cat("\n")
}
