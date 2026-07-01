###supplemental Figure 7: PD data for EYH, LCZ in other organs###
library(tidyverse)
library(janitor)

sup_pd <- read_csv("data/EYH LCZ PD-2024-01( PD-2023-28) report - machine-readable-plates.csv") %>%
  as_tibble() %>%
  clean_names() 

#pivot longer, and setup for dataviz
sup_pd_longer <- sup_pd %>%
  dplyr::filter(compound %in% c("EYH", "LCZ")) %>%
  pivot_longer(cols=c(plasma_conc, spleen_conc, quadricep_conc), names_to = "organ", values_to = "concentration") %>%
  mutate(x_pos = case_when(
    organ == "spleen_conc"    ~ 1,
    organ == "quadricep_conc" ~ 2.5),
    my_colors = if_else(grepl("EYH", compound), eyh_col, lcz_col),
    my_pch=1)

sup_pd_longer$concentration <- as.numeric(unlist(sup_pd_longer$concentration))

#using the MW found in the original PD Google Sheet
eyh_mw = 552.46
lcz_mw = 590.76

#For non-plasma unit conversions, used the molarity formula found in the original Google Sheets: concentration*(0.000000001)/$V$2*1000 -- where $V$2 is the MW
sup_pd_longer <- sup_pd_longer %>%
  mutate(molarity = case_when(
    #when not plasma, since different units)
    #ultimately multiply by 1e6 to go from M to uM
    compound == 'EYH' ~ (concentration * 1e-9 / eyh_mw * 1000) * 1e6,
    compound == 'LCZ' ~ (concentration * 1e-9 / lcz_mw * 1000) * 1e6,
  ))

pd_style <- list(
  xlims = c(0.5, 3), ylims = c(0.1, 115), log = 'y',
  yats = rep(1:9, times=3) * rep(10^(-1:2), each=9),
  ybigs = c(0.1, 1, 10, 100), ybiglabs = c(0.1, 1, 10, 100),
  xvar = "x_pos", yvar = "molarity",
  bartype = 'segment', barwidth = 0.2, boxwidth = NA,
  pchvar = "my_pch", colorvar = "my_colors", jitamt = 0.25,
  xbigs = c(1,2.5), xbiglabs = c("spleen", "quadriceps"),
  ylab = "", mar = c(2,3, 4, 1), yaxcex = 1.1, xaxcex=1.1
)


if (!dir.exists("display_items")) {
  dir.create("display_items")
}

png('display_items/figure_s7.png', width = 6.3, height = 3, units = "in", res = 600)
layout_matrix <- matrix(c(1, 2), nrow = 1)
layout(layout_matrix, widths = c(0.5, 0.5))    

do.call(dviz, c(list(tbl = sup_pd_longer %>% dplyr::filter(compound == "EYH")), pd_style)); abline(h = 1.9, lty = 2)
mtext("compound (\U003BCM)", side = 2, outer = F, line = 2, cex = 1.1, font=1)
mtext("EYH", side = 3, outer = F, line = 2, cex = 1.5, font=1)
do.call(dviz, c(list(tbl = sup_pd_longer %>% dplyr::filter(compound == "LCZ")), pd_style)); abline(h = 0.5, lty = 2)
mtext("LCZ", side = 3, outer = F, line = 2, cex = 1.5, font=1)

dev.off()