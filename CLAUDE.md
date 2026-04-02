# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a scientific research codebase for a manuscript analyzing EYH and LCZ compounds as reducers of prion protein (PrP/PRNP) expression. Data modalities include high-throughput GFP-based screening, flow cytometry, qPCR, TMT proteomics (with phosphoproteomics), ELISA, immunofluorescence microscopy, and in vivo pharmacokinetics.

## Running Scripts

Scripts are standalone R files with no build system. Run individually:

```bash
Rscript src/<script_name>.R
```

Or open in RStudio and run interactively. There is no Makefile or orchestration layer.

**Critical path issue**: Most scripts call `setwd('~/d/sci/src/kd_moa')` at the top, meaning they expect to run from a sibling directory (`kd_moa`) that is not this repo. Data reads like `read_tsv('output/tmt_18h_eyh.tsv')` are relative to that directory. The `output/` and `data/` folders visible in this repo correspond to that working directory.

Scripts also `source('../helper.R')` from `~/d/sci/src/helper.R` — a shared utility file outside this repo providing formatting helpers (percent formatting, confidence intervals, alpha color handling).

## Dependencies

No lockfile exists. Install R packages manually:

```r
install.packages(c("tidyverse", "janitor", "openxlsx", "magick", "drc",
                   "ggridges", "patchwork", "googlesheets4"))

# Bioconductor packages:
BiocManager::install(c("flowCore", "ggcyto", "fgsea", "org.Mm.eg.db",
                       "clusterProfiler", "enrichplot"))

# CRAN:
install.packages("msigdbr")
```

## Code Architecture

### Script purposes
- `figure_2_and_s1.R` — Primary screening figure. Defines a large custom base-R plotting function `dviz()` used throughout for publication-quality plots with KS test statistics. Also processes high-throughput screening data from `data/8pt_*.csv`.
- `figure_3a_eyhlcz.R` — qPCR analysis of Prnp transcript reduction.
- `figure_4_flow_eyh-lcz.R` — Flow cytometry using flowCore/ggcyto for U251MG and HEK293 cells.
- `figure_5.R` — ELISA data for PrP quantification in brain/colon tissue.
- `figure_s2_eyhlcz.R` — TRIM66 peptide-level TMT proteomics analysis.
- `figure_s3_cas9_flow.R` — CRISPR/Cas9 validation flow cytometry across multiple cell lines.
- `figure_S4.R` — Pharmacokinetics (PD) data for EYH and LCZ.
- `gsea_analysis.R` — GSEA on TMT proteomics; reads `output/tmt_18h_{eyh,lcz,y_320}.tsv`; writes GSEA result tables to `output/`.
- `gsea_enrichment_plots.R` — Network and tree plots of GSEA results using clusterProfiler/enrichplot.
- `format_tables_for_docs.R` — Cleans GSEA output for manuscript tables.
- `replace_compound_ids.R` — Maps compound IDs to InChI keys using `output/moa.tsv`.
- `parse_uniprot_mouse.py` — One-off Python script to parse a UniProt XML dump into `output/uniprot_data_mouse.tsv`.

### Data flow
```
data/8pt_*.csv (raw screening)
  → figure_2_and_s1.R → display_items/figure_2.png

data/tmt/*.xlsx (raw TMT proteomics)
  → [prior preprocessing, not in this repo] → output/tmt_18h_{eyh,lcz,y_320}.tsv
  → gsea_analysis.R → output/gsea_*.tsv
  → gsea_enrichment_plots.R → display_items/

output/tmt_18h_*.tsv → figure_s2_eyhlcz.R → display_items/
```

### Color conventions (used consistently across all figures)
- EYH: `#FE7F03` (orange)
- LCZ: `#7B01FC` (purple)
- DMSO (control): `darkgrey`
- Y-320 (positive control): used for comparison
- GFP signal: `#00FF00`; GFP-GPI signal: `#00CDFF`; PrP: `#FF0000`

### Key output columns in TMT tables (`output/tmt_18h_*.tsv`)
- `l2fc`: log2 fold change vs. DMSO
- `p_ebm`: empirical Bayes moderated p-value
- Downstream scripts compute `padj` via `p.adjust(p_ebm, method="fdr")`
