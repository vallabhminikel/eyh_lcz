library(openxlsx)

###figure S2###
dmso_col <- "darkgrey"
eyh_col <- "#FE7F03"
lcz_col <- '#7B01FC'
y320_col <- "#7EF4CC"

#figure S2A####
pep = read.xlsx('data/tmt/1669_P1_peptide.xlsx',startRow=3) %>%
  as_tibble() %>%
  clean_names()

# Trim66 mouse: https://www.uniprot.org/uniprotkb/Q924W6/entry
pep %>%
  filter(gene_symbol=='Trim66') %>%
  mutate(pepno = row_number()) %>%
  select(gene_symbol, pepno, peptide_sequence, dmso_1:y_320_4) %>%
  pivot_longer(cols=dmso_1:y_320_4) %>%
  mutate(name=gsub('y_320','y320',name)) %>%
  mutate(condition = gsub('_.*','',name)) %>%
  mutate(replicate = gsub('.+_','',name)) %>%
  rename(intensity=value) %>%
  select(-name) %>%
  mutate(condx = case_when(condition=='dmso' ~ 1,
                           condition=='eyh' ~ 2,
                           condition=='lcz' ~ 3,
                           condition=='y320' ~ 4)) %>%
  mutate(x = (pepno-1)*5+condx) -> trimpeps

trimpeps %>%
  distinct(pepno, peptide_sequence) %>%
  mutate(midx = (pepno-1)*5 + 2.5) -> whichpeps

trimpeps <- trimpeps %>%
  mutate(
    trim_colors = case_when(
      condition == 'dmso' ~ dmso_col,
      condition == 'eyh'  ~ eyh_col,
      condition == 'lcz'  ~ lcz_col,
      condition == 'y320' ~ y320_col),
    trim_pch = 1
  )

if (!dir.exists("display_items")) {
  dir.create("display_items")
}

resx=600
png('display_items/figure_S2A.png',width=3.5*resx, height=2*resx, res=resx)
par(mar=c(1,2,1,1))

smry= dviz(tbl=trimpeps, xlims=c(0.5, 9.5), ylims=c(0, 60), xvar = 'x', yvar = 'intensity', 
           xcols = c('condition','pepno','peptide_sequence'), yaxcex=0.8, ylabcex=0.8, ylab = 'intensity', ylabline=1.6, colorvar = 'trim_colors',
           log='', xbigs=c(1:4,6:9), xbiglabs=NA, ybigs=0:5*10,xlwds=0, barwidth = 0.33, mar=c(1,3,1,1), pchvar = 'trim_pch')
mtext(side=1, at=smry$x, text=toupper(smry$condition), cex=0.6)
mtext(side=3, line=0, at=whichpeps$midx, text=whichpeps$peptide_sequence, cex=0.45)
dev.off()


#figure S2B####
trim_eyhlcz_qpcr <- read_csv("data/s2_trim_qpcr.csv") %>%
  as_tibble %>%
  clean_names()

trim_eyhlcz_qpcr_long <- trim_eyhlcz_qpcr %>%
  pivot_longer(cols = c(dmso, eyh, lcz), 
               names_to = "treatment", 
               values_to = "value")

trim_eyhlcz_qpcr_long <- trim_eyhlcz_qpcr_long %>%
  mutate(
    my_colors = case_when(
      treatment == 'dmso' ~ dmso_col,
      treatment == 'eyh' ~ eyh_col,
      treatment == 'lcz' ~ lcz_col
    ),
    my_pch = 1,
    x_position = case_when(
      treatment == 'dmso' ~ 0.95,
      treatment == 'eyh' ~ 1.5,
      treatment == 'lcz' ~ 2.05
    ),
    value = value * 100
  )
resx = 600

png('display_items/figure_S2B.png',
    width = 3 * resx,
    height = 2.5 * resx,
    res = resx)

trim_eyhlcz_qpcr_dataviz <- dviz(tbl=trim_eyhlcz_qpcr_long,xvar='x_position', yvar='value',xlims=c(0.6,2.4),ylims=c(0,250),log='',
                                 ybigs=c(0,50,100,150,200,250),ybiglabs=c(0,50,100,150,200,250),bartype='segment',barwidth=0.2,
                                 boxwidth=NA,pchvar='my_pch',colorvar='my_colors',xbigs=c(0.95,1.5,2.05),xbiglabs=c('DMSO','EYH','LCZ'), 
                                 ylab=expression(paste("Normalized", italic("Trim66"), " (%)")),test='t.test',control_group ='dmso',mar=c(2.5,5,3,2),
                                 xaxcex = 0.8)
dev.off()


prnp_eyhlcz_qpcr <- read_csv("data/s2_prnp_qcr.csv") %>%
  as_tibble %>%
  clean_names()

prnp_eyhlcz_qpcr_long <- prnp_eyhlcz_qpcr %>%
  pivot_longer(cols = c(dmso, eyh, lcz), 
               names_to = "treatment", 
               values_to = "value")

prnp_eyhlcz_qpcr_long <- prnp_eyhlcz_qpcr_long %>%
  mutate(
    my_colors = case_when(
      treatment == 'dmso' ~ dmso_col,
      treatment == 'eyh' ~ eyh_col,
      treatment == 'lcz' ~ lcz_col
    ),
    my_pch = 1,
    x_position = case_when(
      treatment == 'dmso' ~ 0.95,
      treatment == 'eyh' ~ 1.5,
      treatment == 'lcz' ~ 2.05
    ),
    value = value * 100
  )

png('display_items/figure_s2b_prnp_qpcr.png',
    width = 3 * resx,
    height = 2.5 * resx,
    res = resx)

prnp_eyhlcz_qpcr_dataviz <- dviz(tbl=prnp_eyhlcz_qpcr_long,xvar='x_position', yvar='value',xlims=c(0.6,2.4),ylims=c(0,250),log='',
                                 ybigs=c(0,50,100,150,200,250),ybiglabs=c(0,50,100,150,200,250),bartype='segment',barwidth=0.2,
                                 boxwidth=NA,pchvar='my_pch',colorvar='my_colors',xbigs=c(0.95,1.5,2.05),xbiglabs=c('DMSO','EYH','LCZ'), 
                                 ylab=expression(paste("Normalized ", italic("Prnp"), " (%)")),test='t.test',control_group ='dmso',mar=c(2.5,5,3,2),
                                 xaxcex = 0.8)
dev.off()


#figure S2C####
trim66_orf_timecourse <- read_csv("data/fig_s2c_orf_transfection.csv") %>%
  as_tibble %>%
  clean_names()

trim66_orf_timecourse_long <- trim66_orf_timecourse %>%
  pivot_longer(
    cols = -time_hours,
    names_to = "treatment_raw",
    values_to = "value"
  ) %>%
  mutate(
    treatment = str_remove_all(treatment_raw, "[_\\.\\d]+$")
  ) %>%
  
  group_by(time_hours, treatment) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  
  select(time_hours, treatment, replicate, value) %>%
  mutate(my_colors = case_when(
    treatment == 'mock' ~ 'black',
    treatment == 'orf' ~ '#1a80bb'
  ),
  my_pch=1)


png('display_items/figure_S2C.png',
    width = 3.5 * resx,
    height = 2.5 * resx,
    res = resx)

trim66_orf_timecourse_dv <- dviz(tbl=trim66_orf_timecourse_long,xvar='time_hours', yvar='value',xlims=c(0,78),ylims=c(0.1,10000),log='y',
                                 ybigs=c(0.1,1,10, 100,1000,10000),ybiglabs=c(expression(10^-1),expression(10^0),expression(10^1), expression(10^2),expression(10^3),expression(10^4)),bartype='segment',barwidth=1,
                                 boxwidth=NA,pchvar='my_pch',colorvar='my_colors',xbigs=c(0, 24, 48, 72),xbiglabs=c(0,24,48,72), 
                                 ylab=expression(paste("Normalized ", italic("Trim66"), " (%)")),mar=c(2.5,5,3,5),
                                 xlab ='hours post-transfection', jitamt=5, legtext = c("Mock", expression(paste(italic("Trim66"), " ORF"))),
                                 legcol = c('black', '#1a80bb'), legpch = c(1, 1),
                                 legcex = 0.7, leglty=0, randseed=6)
dev.off()

#figure S2D####
trim66_si_timecourse <- read_csv("data/fig_s2d_sirna_timecourse.csv") %>%
  as_tibble %>%
  clean_names()


trim66_si_timecourse <- trim66_si_timecourse %>%
  pivot_longer(
    cols = -time_hours,
    names_to = "treatment_raw", 
    values_to = "value"
  ) %>%
  mutate(
    treatment = str_remove_all(treatment_raw, "[_\\.\\d]+$")
  ) %>%
  
  group_by(time_hours, treatment) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  
  select(time_hours, treatment, replicate, value) %>%
  mutate(my_colors = case_when(
    treatment == 'ntc' ~ alpha('black', ci_alpha),
    treatment == 'sirna' ~ alpha('#a00000',ci_alpha)
  ),
  my_pch=1,
  value = value*100
  )

png('display_items/figure_S2D.png',
    width = 3.5 * resx,
    height = 2.5 * resx,
    res = resx)

trim66_orf_timecourse_dv <- dviz(tbl=trim66_si_timecourse,xvar='time_hours', yvar='value',xlims=c(0,78),ylims=c(0,175),log='',
                                 ybigs=c(0,50,100,150),ybiglabs=c(0,50,100,150),bartype='segment',barwidth=1,
                                 boxwidth=NA,pchvar='my_pch',colorvar='my_colors',xbigs=c(0, 24, 48, 72),xbiglabs=c(0,24,48,72), 
                                 ylab=expression(paste("Normalized ", italic("Trim66"), " (%)")),mar=c(2.5,5,3,6),
                                 xlab ='hours post-transfection', jitamt=0.8, legtext = c("Non-targeting \ncontrol", expression(paste(italic("Trim66"), " siRNA"))),
                                 legcol = c('black', '#a00000'), legpch = c(1, 1), legcex = 0.7, leglty=0)

dev.off()