###FIGURE 5####
library(tidyverse)
library(janitor)

eyh_col = "#FE7F03"
lcz_col = "#7B01FC"

read_elisa_data <- function(sheet_id, treat_color) {
  read_csv(sheet_id) %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "treatment_group", values_to = "intensity") %>%
    drop_na(intensity) %>%
    mutate(
      my_colors = if_else(treatment_group == "Vehicle", "darkgrey", treat_color),
      x_pos = if_else(treatment_group == "Vehicle", 1, 1.5),
      my_pch = 1
    )
}
####Figure 5A - PRNP brain ELISA####
#note: each compound has its own vehicle control 
brain_eyh_sheet_id = "data/csELISA_219 analyzed data - brain_elisa_eyh.csv"
brain_lcz_sheet_id = "data/csELISA_219 analyzed data - brain_elisa_lcz.csv"

brain_eyh <- read_elisa_data(brain_eyh_sheet_id, eyh_col)
brain_lcz <- read_elisa_data(brain_lcz_sheet_id, lcz_col)

brain_combined <- bind_rows(
  brain_eyh %>% mutate(
    x_pos = if_else(treatment_group == "Vehicle", 1, 2),
    x_lab = if_else(treatment_group == "Vehicle", "EYH Vehicle", "EYH")
  ),
  brain_lcz %>% mutate(
    x_pos = if_else(treatment_group == "Vehicle", 3, 4),
    x_lab = if_else(treatment_group == "Vehicle", "LCZ Vehicle", "LCZ")
  )
)

#ks test was done previously, but difficult to implement into a graph with 2 different controls. instead, adding * to what was significant in ks test when dviz function done with 1 control.
brain_args <- list(
  xlims = c(0.7, 1.8), ylims = c(0, 130), xvar = "x_pos", yvar = "intensity",
  xbigs = c(1, 1.5), ybigs = seq(0, 130, 25), colorvar = "my_colors",
  barwidth = 0.2, bartype = "bar", ylab = "", log = '', 
  test = NA, control_group = "Vehicle", pchvar = 'my_pch', mar = c(2, 3, 4, 1),
  yaxcex = 1.2, xaxcex=1.2
)

combined_args_ <- brain_args
combined_args_$xlims <- c(0.5, 4.5)
combined_args_$xbigs <- 1:4
combined_args_$xbiglabs <- c("Vehicle", "EYH", "Vehicle", "LCZ")

colon_eyh_sheet_id = "data/6D11_elisa_023_analyzed - colon_elisa_eyh.csv"
colon_lcz_sheet_id = "data/6D11_elisa_023_analyzed - colon_elisa_lcz.csv"

colon_eyh <- read_elisa_data(colon_eyh_sheet_id, eyh_col)
colon_lcz <- read_elisa_data(colon_lcz_sheet_id, lcz_col)

colon_combined <- bind_rows(
  colon_eyh %>% mutate(
    x_pos = if_else(treatment_group == "Vehicle", 1, 2),
    x_lab = if_else(treatment_group == "Vehicle", "Vehicle", "EYH")
  ),
  colon_lcz %>% mutate(
    x_pos = if_else(treatment_group == "Vehicle", 3, 4),
    x_lab = if_else(treatment_group == "Vehicle", "Vehicle", "LCZ")
  )
)

#using the same dataviz arguments as figure 5A, except with change of y-axis ticks
colon_args <- brain_args
colon_args$ylims <- c(0, 200)
colon_args$ybigs <- seq(0, 200, 50)


####Figure 5C - PK, PD data####
fig5c <- read_csv("data/EYH LCZ PD-2024-01( PD-2023-28) report - machine-readable-plates.csv") %>%
  as_tibble() %>%
  clean_names() 

#pivot longer, and setup for dataviz
fig5c_longer <- fig5c %>%
  filter(compound %in% c("EYH", "LCZ")) %>%
  pivot_longer(cols=c(brain_conc, colon_conc, plasma_conc), names_to = "organ", values_to = "concentration") %>%
  mutate(x_pos = case_when(
    organ == "brain_conc"    ~ 1,
    organ == "colon_conc"    ~ 3,
    organ == "plasma_conc"   ~ 5),
    my_colors = if_else(grepl("EYH", compound), eyh_col, lcz_col),
    my_pch=1)

fig5c_longer$concentration <- as.numeric(unlist(fig5c_longer$concentration))

#using the MW found in the original PD Google Sheet
eyh_mw = 552.46
lcz_mw = 590.76

#For unit conversions, used the molarity formula found in the original Google Sheets: concentration*(0.000000001)/$V$2*1000 -- where $V$2 is the MW
fig5c_longer <- fig5c_longer %>%
  mutate(molarity = case_when(
    #when plasma, starting concentration is ng/mL. multiply by 10e3 to convert mL to L, divide by 10e9 to convert ng to g
    #ultimately multiply by 1e6 to go from M to uM
    organ == "plasma_conc" & compound == 'EYH' ~ (concentration * 1e3 / 1e9 / eyh_mw) * 1e6,
    organ == "plasma_conc" & compound == 'LCZ' ~ (concentration * 1e3 / 1e9 / lcz_mw) * 1e6,
    #ultimately multiply by 1e6 to go from M to uM
    organ != "plasma_conc" & compound == 'EYH' ~ (concentration * 1e-9 / eyh_mw * 1e3) * 1e6,
    organ != "plasma_conc" & compound == 'LCZ' ~ (concentration * 1e-9 / lcz_mw * 1e3) * 1e6
  ))

fig5c_eyh <- fig5c_longer %>% filter(compound=="EYH")
fig5c_lcz <- fig5c_longer %>% filter(compound=="LCZ")

fig5c_style <- list(
  xlims = c(0.5, 5.5), ylims = c(0.1, 100), log = 'y',
  yats = rep(1:9, times=3) * rep(10^(-1:2), each=9),
  ybigs = c(0.1, 1, 10, 100), ybiglabs = c(0.1, 1, 10, 100),
  xvar = "x_pos", yvar = "molarity",
  bartype = 'segment', barwidth = 0.2, boxwidth = NA,
  pchvar = "my_pch", colorvar = "my_colors", jitamt = 0.25,
  xbigs = c(1, 3, 5), xbiglabs = c("brain", "colon","plasma"),
  ylab = "", mar = c(2,3, 4, 1),
  yaxcex = 1.2, xaxcex=1.2
)


if (!dir.exists("display_items")) {
  dir.create("display_items")
}

###combine into single display item####
alphabet = 2.5
axis_size = 1

#aiming for 4 panels total (2 ELISA, 2 PK)
layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 1)
png('display_items/figure_5.png', width = 8, height = 3, units = "in", res = 600)
layout(layout_matrix, widths = c(1, 1, 0.75, 0.75))
par(oma = c(3, 4, 3, 1), ps = 10)

#figure 5A: Brain
do.call(dviz, c(list(tbl = brain_combined), combined_args_))
mtext("A", side = 3, adj = -0.2, line = 2, font = 1, cex = 1.5)
mtext("Brain", side = 3, line = 1, cex = 1.2)
mtext("Normalized PrP (%)", side = 2, line = 2.5)

#figure 5B: Colon
colon_args_final <- combined_args_
colon_args_final$ylims <- c(0, 200)
colon_args_final$ybigs <- seq(0, 200, 50)

do.call(dviz, c(list(tbl = colon_combined), colon_args_final))
text(x = 2, y = 195, labels = "*", cex = 1.5, font = 1)
mtext("B", side = 3, adj = -0.2, line = 2, font = 1, cex = 1.5)
mtext("Colon", side = 3, line = 1, cex = 1.2)
mtext("Normalized PrP (%)", side = 2, line = 2.5)

#figure 5C: PK (making these as separate small panels)
do.call(dviz, c(list(tbl = fig5c_eyh), fig5c_style)); abline(h = 1.9, lty = 2)
mtext("C", side = 3, adj = -0.4, line = 2, font = 1, cex = 1.5)
mtext("EYH", side = 3, line = 1,cex=1.2)
mtext("compound (\U003BCM)", side = 2, outer = F, line = 2, cex = 1, font=1)
do.call(dviz, c(list(tbl = fig5c_lcz), fig5c_style)); abline(h = 0.5, lty = 2)
mtext("LCZ", side = 3, line = 1, cex=1.2)

dev.off()