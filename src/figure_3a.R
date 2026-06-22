####Figure 3####
dmso_col <- "#777777"
eyh_col <- "#FE7F03"
lcz_col <- '#7B01FC'

#plot qPCR that Prnp-lowering effects of EYH/LCZ are not transcriptional
fig3 <- read_csv("data/qpcr-3-10-3-11.csv") %>%
  as_tibble %>%
  clean_names() %>%
  select(-average_dmso)

#all data are given a number for biological replicate. want to plot all of them, so don't care which number they are. 
#just need to strip their numerical ID using gsub
fig3 <- fig3 %>%
  mutate(
    treatment_group = gsub(" [0-9]", "", tx),
    my_colors = case_when(
      treatment_group == "DMSO" ~ dmso_col,
      treatment_group == "EYH"  ~ eyh_col,
      treatment_group == "LCZ"  ~ lcz_col
    ),
    x_position = case_when(
      treatment_group == "DMSO" ~ 1,
      treatment_group == "EYH"  ~ 1.5,
      treatment_group == "LCZ"  ~ 2
    ),
    my_pch = 21
  ) %>%
  filter(treatment_group != 'Y-320') %>%
  rename(x = x_position, intensity = normalized_to_dmso) %>%
  group_by(date, treatment_group, tx, x, my_colors, my_pch) %>%
  summarize(.groups='keep',
            mean_intensity = mean(intensity)) %>%
  ungroup()


resx = 600

png('display_items/figure_3A.png',
    width = 3.5 * resx,
    height = 2 * resx,
    res = resx)

#use KS test
result = dviz(tbl=fig3, xlims=c(0.7,2.3), ylims=c(0,130), xvar="x", yvar="mean_intensity", 
              xbigs=c(1,1.5,2), xbiglabs=c("DMSO","EYH","LCZ"), ybigs=seq(0,130,25), 
              xcols = c('treatment_group'),
              colorvar="my_colors", barwidth=0.075, bartype="bar", 
              ylab=expression(paste("Normalized ", italic("Prnp"), " (%)")), 
              log='', test='ks.test', control_group="DMSO", pchvar="my_pch", mar=c(1.5,4,1,1))

abline(h = 100, lty = 2, lwd = 1.5) 
dev.off()