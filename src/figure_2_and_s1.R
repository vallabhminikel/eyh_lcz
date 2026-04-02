library(tidyverse)
library(janitor)
library(openxlsx)
library(magick)
library(drc); select=dplyr::select; summarize=dplyr::summarize
setwd('~/d/sci/src/kd_moa')
source('../helper.R')


gpicol = '#00CDFF'
comcol = '#008888'

defcol = '#A9A9A9'
prpcol = '#FF0000'
gfpcol = '#00FF00'
toxcol = '#0000FF'



dviz = function(tbl,
                xlims,
                ylims,
                xvar,
                yvar,
                colorvar=NULL,
                pchvar=NULL,
                pchbg='#FFFFFF',
                pchout=NULL,
                pchcex=1,
                xcols=character(0), # columns that group with x TBD
                xats=xbigs,
                xbigs,
                xlwds=1,
                xbiglabs=xbigs,
                xaxcex=1,
                xlabcex=1,
                xlabline=1.5,
                yats=ybigs,
                ybigs,
                ylwds=1,
                ybiglabs=ybigs,
                yaxcex=1,
                ylabcex=1,
                ylabline=2.25,
                log,
                mar=c(3,3,1,1),
                jitamt=0.1,
                randseed=1,
                boxwidth=NA,
                barwidth=NA,
                bartype='segment',
                polygon=NA,
                arrowlength=0.05,
                test=NA,
                control_group=NA,
                xlab='',
                ylab='',
                legtext=NULL,
                legcol='#000000',
                legtextcol=legcol,
                leglty=1,
                legpch=20,
                leglwd=1,
                legcex=1,
                legloc=NULL,
                crosshairs=F
) {
  
  if (is.null(pchvar)) {
    pchvar='pch'
    tbl$pch = 19
  }
  
  if (is.null(colorvar)) {
    colorvar='color'
    tbl$color = '#000000'
  }
  
  tbl %>%
    mutate(x=!!as.name(xvar), y=!!as.name(yvar), color=!!as.name(colorvar), pch=!!as.name(pchvar)) %>%
    select(x, y, color, pch, all_of(xcols)) -> tbl
  
  if (!crosshairs) {
    xcols = c('x', xcols)
  }
  
  tbl %>%
    group_by(color, across(all_of(xcols))) %>%
    summarize(.groups='keep',
              n = n(),
              mean = mean(y),
              l95 = mean(y) - 1.96 * sd(y) / sqrt(n()),
              u95 = mean(y) + 1.96 * sd(y) / sqrt(n()),
              median = median(y),
              q25 = quantile(y, .25),
              q75 = quantile(y, .75),
              cv = sd(y) / mean(y),
              x_mean = mean(x),
              x_l95  = mean(x) - 1.96 * sd(x) / sqrt(n()),
              x_u95  = mean(x) + 1.96 * sd(x) / sqrt(n())) %>%
    ungroup() %>%
    mutate(l95 = ifelse(l95 < 0 & log %in% c('xy','y'), min(ylims), l95)) %>%
    #mutate(x = case_when(crosshairs ~ x_mean, 
    #                     TRUE ~ x)) %>%
    arrange(x) -> tbl_smry
  
  if (crosshairs) {
    tbl_smry$x = tbl_smry$x_mean
  }
  
  par(mar=mar)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log=log)
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  if (!is.null(xats)) {
    axis(side=1, at=xats, tck=-0.025, lwd.ticks=xlwds, labels=NA)
  }
  if (!is.null(xbigs)) {
    axis(side=1, at=xbigs, tck=-0.05, lwd.ticks=xlwds, labels=NA)
    axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbiglabs, cex.axis=xaxcex)
  }
  mtext(side=1, line=xlabline, text=xlab, cex=xlabcex)
  if (!is.null(yats)) {
    axis(side=2, at=yats, tck=-0.025, labels=NA)
  }
  if (!is.null(ybigs)) {
    axis(side=2, at=ybigs, tck=-0.05, labels=NA)
    axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.3, labels=ybiglabs, cex.axis=yaxcex)
  }
  mtext(side=2, line=ylabline, text=ylab, cex=ylabcex)
  if (crosshairs) {
    jitamt = 0
  }
  
  if (crosshairs) {
    segments(x0=tbl_smry$x_l95, x1=tbl_smry$x_u95, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    segments(x0=tbl_smry$x_mean, y0=tbl_smry$l95, y1=tbl_smry$u95, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(barwidth)) {
    if (bartype=='segment') {
      segments(x0=tbl_smry$x-barwidth, x1=tbl_smry$x+barwidth, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    } else if (bartype=='bar') {
      rect(xleft=tbl_smry$x-barwidth, xright=tbl_smry$x+barwidth, ybottom=rep(0,nrow(tbl_smry)), ytop=tbl_smry$mean, col=tbl_smry$color, border=NA)
    }
    arrows(x0=tbl_smry$x, y0=tbl_smry$l95, y1=tbl_smry$u95, code=3, angle=90, length=arrowlength, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(boxwidth)) {
    rect(xleft=tbl_smry$x-boxwidth, xright=tbl_smry$x+boxwidth, ybottom=tbl_smry$q25, ytop=tbl_smry$q75, border=tbl_smry$color, lwd=1.5, col=NA)
    segments(x0=tbl_smry$x-boxwidth, x1=tbl_smry$x+boxwidth, y0=tbl_smry$median, col=tbl_smry$color, lwd=1.5)
  }
  
  set.seed(randseed)
  if (!is.null(pchout)) {
    points(x=jitter(tbl$x,amount=jitamt), y=tbl$y, col=pchout, pch=tbl$pch, bg=pchbg)
  } else {
    points(x=jitter(tbl$x,amount=jitamt), y=tbl$y, col=alpha(tbl$color,ci_alpha), pch=tbl$pch, bg=pchbg) 
  }
  
  
  
  if (!is.na(polygon)) {
    for (clr in unique(tbl_smry$color)) {
      if (polygon=='iqr') {
        subs = subset(tbl_smry, color==clr & !is.na(q25) & !is.na(q75))
        points( x=  subs$x, y=  subs$median, type='l', lwd=2, col=subs$color)
        polygon(x=c(subs$x, rev(subs$x)), y=c(subs$q25, rev(subs$q75)), col=alpha(subs$color, ci_alpha), border=NA)
      } else if (polygon=='ci') {
        subs = subset(tbl_smry, color==clr & !is.na(l95) & !is.na(u95))
        points( x=  subs$x, y=  subs$mean, type='l', lwd=2, col=subs$color)
        polygon(x=c(subs$x, rev(subs$x)), y=c(subs$l95, rev(subs$u95)), col=alpha(subs$color, ci_alpha), border=NA)
      }
    }
  }
  
  if (!is.na(test)) {
    testfun = get(test) # e.g. ks.test
    control_color = control_group
    tbl_smry$p = as.numeric(NA)
    for (i in 1:nrow(tbl_smry)) {
      this_x = tbl_smry$x[i]
      this_rows = tbl$x == this_x & tbl$color == tbl_smry$color[i]
      ctrl_rows = tbl$x == this_x & tbl$color == control_group
      test_obj = suppressWarnings(testfun(tbl$y[this_rows],  tbl$y[ctrl_rows]))
      tbl_smry$p[i] = test_obj$p.value
    }
    tbl_smry$p_symb = ''
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.05] = '*'
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.01] = '**'
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.001] = '***'
    text(x=tbl_smry$x[tbl_smry$color != control_color], y=max(ylims)*.95, labels=tbl_smry$p_symb[tbl_smry$color != control_color])
  }
  
  if (!is.null(legtext)) {
    if (is.null(legloc)) {
      par(xpd=T)
      legend(x=max(xlims),y=max(ylims),legtext,col=legcol,text.col=legtextcol,pch=legpch,lwd=leglwd, cex=legcex,lty=leglty, bty='n')
      par(xpd=F)
    } else {
      if (length(legloc)==1) {
        legend(legloc,legtext,col=legcol,text.col=legtextcol,pch=legpch,lwd=leglwd, cex=legcex,lty=leglty, bty='n')
      } else if (length(legloc)==2) {
        legend(x=legloc[1],y=legloc[2],legtext,col=legcol,text.col=legtextcol,pch=legpch,lwd=leglwd, cex=legcex,lty=leglty, bty='n')
      }
    }
  }
  
  return(tbl_smry)
  
}

# these files only contain one row per sample_name + wells_concentration combination
# the screen was in duplicate. presumably these are already averaged between duplicates.
ht_gfp = read_csv('data/hit_table_gfp.csv', col_types=cols()) %>%
  clean_names() %>%
  mutate(cell_line = 'N2a-GFP')
ht_gpi = read_csv('data/hit_table_gfp-gpi.csv', col_types=cols()) %>%
  clean_names() %>%
  mutate(cell_line = 'N2a-GFP-GPI')

# what about this? from here: https://drive.google.com/drive/folders/1s-QXkSKnESEA2yfV00LQe80nKIuD3V1W?usp=drive_link
dat = read_csv('data/gfp_well-level_correctedScoresOutput_data.csv') %>%
  clean_names()
dat = read_csv('data/gfp-gpi_well-level_correctedScoresOutput_data.csv') %>%
  clean_names()
dat %>%
  group_by(inchi_key, concentration) %>%
  summarize(.groups='keep', n=n())
# good - this has the dups
# it appears that well_type == 'R1' is GFP siRNA and well_type == 'R2' is PrP siRNA

dat %>%
  mutate(condition = case_when(well_type=='R2' ~ 'poscon',
                               well_type=='NC' ~ 'negcon')) %>%
  filter(condition %in% c('poscon','negcon')) %>%
  mutate(readout = median_membrane_intensity_integrated_intensity_pr_p_value) %>%
  group_by(condition) %>%
  summarize(.groups='keep',
            mean_readout = mean(readout, na.rm=T),
            sd_readout = sd(readout, na.rm=T)) -> zpt

z_prime = 1 - (3*zpt$sd_readout[zpt$condition=='poscon'] + 3*zpt$sd_readout[zpt$condition=='negcon']) / 
  abs(zpt$mean_readout[zpt$condition=='poscon'] - zpt$mean_readout[zpt$condition=='negcon'])


dat %>%
  mutate(condition = case_when(well_type=='R2' ~ 'poscon',
                               well_type=='NC' ~ 'negcon')) %>%
  filter(condition %in% c('poscon','negcon')) %>%
  mutate(readout = median_cytoplasm_intensity_mean_intensity_pr_p_value) %>%
  group_by(condition) %>%
  summarize(.groups='keep',
            mean_readout = mean(readout, na.rm=T),
            sd_readout = sd(readout, na.rm=T)) -> zpt

z_prime = 1 - (3*zpt$sd_readout[zpt$condition=='poscon'] + 3*zpt$sd_readout[zpt$condition=='negcon']) / 
  abs(zpt$mean_readout[zpt$condition=='poscon'] - zpt$mean_readout[zpt$condition=='negcon'])


dat %>%
  mutate(condition = case_when(well_type=='R2' ~ 'poscon',
                               well_type=='NC' ~ 'negcon')) %>%
  filter(condition %in% c('poscon','negcon')) %>%
  select(condition, matches('pr_p')) %>%
  pivot_longer(cols=matches('pr_p')) %>%
  group_by(condition, name) %>%
  summarize(.groups='keep',
            mean_readout = mean(value, na.rm=T),
            sd_readout = sd(value, na.rm=T)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(.groups='keep',
            zprime = 1 - (3*sd_readout[condition=='poscon'] + 3*sd_readout[condition=='negcon']) / 
              abs(mean_readout[condition=='poscon'] - mean_readout[condition=='negcon'])) %>%
  ungroup() -> r2_zprimes

write_tsv(r2_zprimes,'display_items/r2_zprimes.tsv')


dat %>%
  mutate(condition = case_when(well_type=='R1' ~ 'poscon',
                               well_type=='NC' ~ 'negcon')) %>%
  filter(condition %in% c('poscon','negcon')) %>%
  select(condition, matches('gfp')) %>%
  pivot_longer(cols=matches('gfp')) %>%
  group_by(condition, name) %>%
  summarize(.groups='keep',
            mean_readout = mean(value, na.rm=T),
            sd_readout = sd(value, na.rm=T)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(.groups='keep',
            zprime = 1 - (3*sd_readout[condition=='poscon'] + 3*sd_readout[condition=='negcon']) / 
              abs(mean_readout[condition=='poscon'] - mean_readout[condition=='negcon'])) %>%
  ungroup() -> r1_zprimes

write_tsv(r1_zprimes,'display_items/r1_zprimes.tsv')



resx=600
png('display_items/replicate_agreement.png', width=6.5*resx, height=3.5*resx, res=resx)

layout_matrix = matrix(1:2, byrow=T, nrow=1)
layout(layout_matrix)
par(mar=c(4,4,1,1))

dat %>%
  filter(!is.na(dat$inchi_key)) %>%
  select(inchi_key, concentration,
         median_cytoplasm_intensity_integrated_intensity_gfp_value,
         median_membrane_intensity_integrated_intensity_pr_p_value) %>%
  mutate(id = paste0(inchi_key, ' ', concentration)) %>%
  mutate(replicate = case_when(duplicated(id) ~ 2,
                                          TRUE ~ 1)) -> readouts_with_replicates
  
readouts_with_replicates %>% 
  select(inchi_key, concentration, replicate, readout=median_cytoplasm_intensity_integrated_intensity_gfp_value) %>%
  pivot_wider(names_from=replicate, values_from=readout) %>%
  rename(r1 = `1`, r2 = `2`) -> replicates

horiz_col = '#776699'
lims = range(c(replicates$r1, replicates$r2))*1.1
plot(replicates$r1, replicates$r2, xaxs='i', yaxs='i', xlim=lims, ylim=lims, axes=F, ann=F, pch=20, col=alpha(gfpcol, ci_alpha))
axis(side=1, at=-7:7*10, lwd=0, cex.axis=0.5)
axis(side=2, at=-7:7*10, las=2, lwd=0, cex.axis=0.5)
axis(side=1, at=-7:7*10, labels = NA, tck=-0.02, pos=0)
axis(side=2, at=-7:7*10, labels = NA, tck=-0.02, pos=0)
mtext(side=1, line=2.2, text='replicate 1')
mtext(side=2, line=2.2, text='replicate 2')
mtext(side=3, line=0, text='cytoplasmic GFP')
abline(a=0, b=1, col=horiz_col)
pearsons = cor.test(replicates$r1, replicates$r2)
rho = pearsons$estimate
p = pearsons$p.value
text(x=min(lims),y=max(lims)*.9,pos=4,cex=0.8,labels=paste0('r = ',formatC(rho, format='f', digits=3)))


readouts_with_replicates %>% 
  select(inchi_key, concentration, replicate, readout=median_membrane_intensity_integrated_intensity_pr_p_value) %>%
  pivot_wider(names_from=replicate, values_from=readout) %>%
  rename(r1 = `1`, r2 = `2`) -> replicates

horiz_col = '#776699'
lims = range(c(replicates$r1, replicates$r2))*1.1
plot(replicates$r1, replicates$r2, xaxs='i', yaxs='i', xlim=lims, ylim=lims, axes=F, ann=F, pch=20, col=alpha(prpcol, ci_alpha))
axis(side=1, at=-7:7*10, lwd=0, cex.axis=0.5)
axis(side=2, at=-7:7*10, las=2, lwd=0, cex.axis=0.5)
axis(side=1, at=-7:7*10, labels = NA, tck=-0.02, pos=0)
axis(side=2, at=-7:7*10, labels = NA, tck=-0.02, pos=0)
mtext(side=1, line=2.2, text='replicate 1')
mtext(side=2, line=2.2, text='replicate 2')
mtext(side=3, line=0, text='membrane PrP')
abline(a=0, b=1, col=horiz_col)
pearsons = cor.test(replicates$r1, replicates$r2)
rho = pearsons$estimate
p = pearsons$p.value
text(x=min(lims),y=max(lims)*.9,pos=4,cex=0.8,labels=paste0('r = ',formatC(rho, format='f', digits=3)))

dev.off()




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

resx=600
png('display_items/trim66_peptides.png',width=3.5*resx, height=2*resx, res=resx)
par(mar=c(1,2,1,1))
smry= dviz(tbl=trimpeps, xlims=c(0.5, 9.5), ylims=c(0, 60), xvar = 'x', yvar = 'intensity', 
     xcols = c('condition','pepno','peptide_sequence'), yaxcex=0.8, ylabcex=0.8, ylab = 'intensity', ylabline=1.6,
     log='', xbigs=c(1:4,6:9), xbiglabs=NA, ybigs=0:5*10,xlwds=0, barwidth = 0.33, mar=c(1,3,1,1))
mtext(side=1, at=smry$x, text=toupper(smry$condition), cex=0.6)
mtext(side=3, line=0, at=whichpeps$midx, text=whichpeps$peptide_sequence, cex=0.45)
dev.off()





tibble(path = list.files('rawimgs/2022_04_07', recursive = T, full.names = T)) %>%
  mutate(tx = gsub('/.*','',gsub('rawimgs/2022_04_07/GFP/', '', path))) %>%
  mutate(tx = case_when(tx == 'PrP' ~ 'PrPsiRNA',
                        tx=='DMSO' ~ 'DMSO',
                        TRUE ~ substr(tx,5,10))) %>%
  mutate(field=gsub('A0[1-4]Z0[1-4]C0[1-4].tif','',gsub('.*/','',path))) %>%
  mutate(channel = gsub('.tif','',str_extract(path, 'C0[1-4].tif'))) %>%
  mutate(stain = case_when(channel == 'C01' ~ 'Hoechst',
                           channel == 'C02' ~ 'Cellmask',
                           channel == 'C03' ~ 'PrP',
                           channel == 'C04' ~ 'GFP')) -> imgs_all

representative_images = tribble(
  ~cpd, ~tx, ~field,
  'PrP_siRNA', 'PrPsiRNA', 'VH97007533_F12.d3_T0001F002L01',
  'DMSO', 'DMSO', 'VH97007533_H12.a2_T0001F002L01',
  'Y-320', 'CCW009', 'VH97007533_G07.a1_T0001F002L01',
  'EYH', 'CLV652', 'VH97007533_G01.a3_T0001F001L01',
  'LCZ', 'LCZ960', 'VH97007527_F01.c4_T0001F002L01'
)

imgs_all %>%
  inner_join(representative_images, by=c('tx','field')) -> imgs


for (field in unique(imgs$field)) {
  
  subs = imgs[imgs$field==field,]
  cat(file=stderr(), paste0('\rNow merging ',subs$tx[1],' field ',field,'...'))
  red = image_read(subs$path[subs$stain=='PrP'])
  grn = image_read(subs$path[subs$stain=='GFP'])
  blu = image_read(subs$path[subs$stain=='Hoechst'])
  red = image_modulate(red, brightness=2500, saturation=100)
  grn = image_modulate(grn, brightness=500, saturation=100)
  blu = image_modulate(blu, brightness=500, saturation=100)
  rgb = image_combine(c(red, grn, blu), colorspace='RGB')
  rgb = image_scale(rgb, '600x600')
  new_path = paste0('display_items/if_example_',subs$cpd[1],'.png')
  image_write(rgb, new_path, format='png')
}


dr8 = read_tsv('output/dr8.tsv')

resx=600
Cairo::CairoPNG('display_items/dose_response.png', width=8.5*resx, height=5*resx, res=resx)
#cairo_pdf('display_items/dose_response.pdf', width=6.5, height=5)
layout_matrix = matrix(1:6, nrow=2, byrow=T)
layout(layout_matrix)

panel = 1

par(mar=c(3,4,3,1))

xlims = c(.003, 30)
ylims = c(0,1.9)


for (these_cells in c('GFP','GFP-GPI')) {


dr8 %>%
  filter(cells==these_cells) %>%
  mutate(tx = case_when(is.na(inchi_key) ~ 'DMSO',
                        inchi_key == 'BWZNJVZTAWBIFG-UHFFFAOYSA-N' ~ 'Y-320',
                        inchi_key == 'LAYAEBGGDTWXAZ-RUZDIDTESA-N' ~ 'LCZ',
                        inchi_key == 'XBJHBEQDHSREQB-QHCPKHFHSA-N' ~ 'EYH',
                        TRUE ~ NA)) %>%
  filter(!is.na(tx) & tx != 'DMSO') %>%
  group_by(tx) %>%
  mutate(nuc_rel = nuclei / mean(nuclei[conc < 0.1]),
         prp_rel = (prp - min(prp)) / (mean(prp[conc < 0.1]) - min(prp)),
         gfp_rel = (gfp - min(gfp)) / (mean(gfp[conc < 0.1]) - min(gfp))) %>%
  ungroup() %>%
  arrange(tx) -> dr8_hits


for (cpd in sort(unique(dr8_hits$tx))) {

  subs = dr8_hits[dr8_hits$tx==cpd,]
  
  nucfit = drm(nuc_rel ~ conc, data=subs, fct=LL.4(fixed=c(b=NA,c=NA,d=1,e=NA)))
  prpfit = drm(prp_rel ~ conc, data=subs, fct=LL.4(fixed=c(b=NA,c=NA,d=1,e=NA)))
  gfpfit = drm(gfp_rel ~ conc, data=subs, fct=LL.4(fixed=c(b=NA,c=NA,d=1,e=NA)))
  
  x = seq(.003,30,.03) # for model fits
  xats = rep(c(1:9),5) * rep(10^(-3:1),each=9)
  xbigs = 10^(-2:1)
  
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=F, axes=F, log='x')
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  mtext(side=1, line=0.5, at=xbigs, text=xbigs, cex=0.6)
  mtext(side=1, line=1.5, text='concentration (µM)', cex=0.8)
  axis(side=2, at=0:20/10, tck=-0.025, labels=NA)
  axis(side=2, at=0:4/2, tck=-0.05, labels=NA)
  mtext(side=2, line=0.75, at=0:3/2, text=percent(0:3/2), las=2, cex=0.6)
  mtext(side=2, line=3, text='residual signal', cex=0.8)
  points(x=subs$conc, y=subs$prp_rel, pch=20, col=alpha(prpcol, ci_alpha))
  y = suppressWarnings(predict(prpfit, newdata=data.frame(concentration=x)))
  points(x=x, y=y, lwd=1, lty=1, type='l', col=prpcol)
  points(x=subs$conc, y=subs$gfp_rel, pch=20, col=alpha(gfpcol, ci_alpha))
  y = suppressWarnings(predict(gfpfit, newdata=data.frame(concentration=x)))
  points(x=x, y=y, lwd=1, lty=1, type='l', col=gfpcol)
  points(x=subs$conc, y=subs$nuc_rel, pch=20, col=alpha(toxcol, ci_alpha))
  y = suppressWarnings(predict(nucfit, newdata=data.frame(concentration=x)))
  points(x=x, y=y, lwd=1, lty=1, type='l', col=toxcol)
  
  mtext(side=3, line=0, text=cpd, cex=0.8)
  
  
  prp_ec50_value = as.numeric(prpfit$coefficients['e:(Intercept)'])
  prp_ec50 = formatC(prp_ec50_value,format='f',digits=2)
  prp_ec50[prp_ec50_value > 20] = '>20' # if >20 set to >20
  
  gfp_ec50_value = as.numeric(gfpfit$coefficients['e:(Intercept)'])
  gfp_ec50 = formatC(gfp_ec50_value,format='f',digits=2)
  gfp_ec50[gfp_ec50_value > 20] = '>20' # if >20 set to >20
  
  cytotox_ec50_value = as.numeric(nucfit$coefficients['e:(Intercept)'])
  cytotox_ec50 = formatC(cytotox_ec50_value,format='f',digits=2)
  cytotox_ec50[ as.numeric(nucfit$coefficients['c:(Intercept)']) > 1] = '>20' # if nuclei count actually increases with dose (c > 1), then no EC50 found
  
  leg = tribble(
    ~color, ~disp,
    prpcol, paste0('PrP EC\U2085\U2080 (\U00B5M): ',prp_ec50),
    toxcol, paste0('Cytotox EC\U2085\U2080 (\U00B5M): ',cytotox_ec50),
    gfpcol, paste0('GFP EC\U2085\U2080 (\U00B5M): ',gfp_ec50)
  )
  
  legend('topleft', leg$disp, col=leg$color, pch=15, cex=0.7, bty='n')
  
  
  panel = panel + 1
}
  
}
dev.off() 






eyh_allcano = read_tsv('output/tmt_18h_eyh.tsv')
lcz_allcano = read_tsv('output/tmt_18h_lcz.tsv')
y320_allcano = read_tsv('output/tmt_18h_y_320.tsv')

rbind(eyh_allcano %>% mutate(tx='EYH'), 
      lcz_allcano %>% mutate(tx='LCZ'),
      y320_allcano %>% mutate(tx='Y-320')) -> allcanos

ylims = c(0, max(-log10(allcanos$p_ebm[allcanos$tx %in% c('EYH','LCZ')]))*1.05)
clipx = 3
xlims = c(-clipx, clipx)

resx=600
png('display_items/figure_2.png', width=6.5*resx, height=7.0*resx, res=resx)

layout_matrix = matrix(c(1:4,
                       rep(5,4),
                       rep(6,4)), nrow=3, byrow=T)
layout(layout_matrix, widths=c(0.25, 1, 1, 1),
       heights = c(1,.5,1))

panel = 1

par(mar=c(3,0,3,0))

plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=4, at=0:20, tck=0.05, labels=NA)
axis(side=4, at=0:20, tck=0.025, lwd=0, las=2, line=-2)
mtext(side=4, line=-3, text='-log10(P value)', cex=0.8)

par(mar=c(3,0,3,1))

for (this_tx in c('EYH','LCZ','Y-320')) {
  allcano = allcanos %>% filter(tx==this_tx)
  
  plot(x=clipdist(allcano$l2fc,-clipx, clipx), y=clipdist(-log10(allcano$p_ebm), 0, max(ylims)), pch=20, axes=F, ann=F, xaxs='i', yaxs='i',
       col=allcano$color,
       xlim=xlims, ylim=ylims)
  axis(side=1, at=c(-clipx, clipx), lwd.ticks=0, labels=NA)
  axis(side=1, at=-clipx:clipx, labels=NA)
  axis(side=1, at=-clipx:clipx, cex.axis=0.8, lwd=0, line=-0.5)
  mtext(side=1, line=1.8, text='log2 fold change', cex=0.8)
  corrected_threshold = 0.05 / nrow(allcano)
  abline(h=-log10(c(0.05, corrected_threshold)), lwd=0.25)
  text(x=c(-3,-3), y=-log10(c(0.05, corrected_threshold))-0.18, pos=4, labels=c('nominal','Bonferroni'), font=3, cex=0.6)
  if (this_tx %in% c('EYH','LCZ')) {
    callouts = -log10(allcano$p_ebm) > 7  & abs(allcano$l2fc > 0.76) | 
      rank(allcano$p_ebm) <= 5 | 
      (this_tx == 'EYH' & -log10(allcano$p_ebm) > 5.32 & abs(allcano$l2fc) > 0.96) |
      (this_tx == 'LCZ' & -log10(allcano$p_ebm) > 5.32 & abs(allcano$l2fc) > 1.00)
  } else {
    callouts = allcano$gene %in% 'Prnp'
  }
  par(xpd=T)
  text(x=clipdist(allcano$l2fc[callouts], -clipx, clipx), y=-log10(allcano$p_ebm[callouts]), 
       pos=ifelse(allcano$l2fc[callouts] < 0, 2, 4), 
       labels=allcano$gene[callouts], font=3, cex=0.5, family='mono')
  par(xpd=F)
  points(allcano$l2fc[allcano$gene=='Prnp'], -log10(allcano$p_ebm[allcano$gene=='Prnp']), pch=1, cex=1.5)
  mtext(side=3, line=0.25, text=this_tx)
  mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1
}

# panel D - localization analysis

locs_raw = read_tsv('output/uniprot_data_mouse.tsv', col_types=cols())
locs_raw %>%
  select(gene, both) %>%
  separate_rows(both, sep=',') %>%
  filter(!is.na(both)) -> locs

locs_to_investigate = c('Cytoplasm','Nucleus','Cell membrane','Secreted')

results = tibble(loc = locs_to_investigate,
                 eyh_es = as.numeric(NA),
                 eyh_p = as.numeric(NA),
                 lcz_es = as.numeric(NA),
                 lcz_p = as.numeric(NA))
for (i in 1:nrow(results)) {
  genes_with_this_loc = locs %>% filter(both==results$loc[i]) %>% pull(gene)
  eyh_allcano %>%
    mutate(bonf = p_ebm < 0.05 / nrow(eyh_allcano)) %>%
    mutate(located = gene %in% genes_with_this_loc) -> fisher_data
  fisher_obj = fisher.test(table(fisher_data[,c('bonf','located')]), alternative='two.sided')
  results$eyh_es[i] = fisher_obj$estimate
  results$eyh_l95[i] = fisher_obj$conf.int[1]
  results$eyh_u95[i] = fisher_obj$conf.int[2]
  results$eyh_p[i] = fisher_obj$p.value
  
  lcz_allcano %>%
    mutate(bonf = p_ebm < 0.05 / nrow(lcz_allcano)) %>%
    mutate(located = gene %in% genes_with_this_loc) -> fisher_data
  fisher_obj = fisher.test(table(fisher_data[,c('bonf','located')]), alternative='two.sided')
  results$lcz_es[i] = fisher_obj$estimate
  results$lcz_l95[i] = fisher_obj$conf.int[1]
  results$lcz_u95[i] = fisher_obj$conf.int[2]
  results$lcz_p[i] = fisher_obj$p.value
}                 
n_tests = 2*nrow(results)
results$eyh_bonf = pmin(1,results$eyh_p * n_tests)
results$lcz_bonf = pmin(1,results$lcz_p * n_tests)


# what else is GPI-anchored
allcanos %>%
  filter(tx %in% c('EYH','LCZ')) %>%
  filter(p_ebm < 0.05 / nrow(eyh_allcano)) %>%
  filter(gene %in% locs$gene[locs$both=='GPI-anchor']) -> gpi_lookup

# panels D & E - vibe coded with Gemini Pro, reviewed by Eric

# Set graphical parameters: adjust left margin for location labels
par(mar=c(3, 14, 3, 2))

# Define colors and plot settings
color_eyh <- "#FE7F03"
color_lcz <- "#7B01FC"
n_locs <- nrow(results)
y_pos <- seq_len(n_locs)
offset <- 0.15  # Stagger groups vertically

# 1. Create a blank plot area with log x-axis
plot(NA, 
     log = "x",                  # Set x-axis to log scale
     xlim = c(0.03, 10),         # Set x-axis limits
     ylim = c(0.5, n_locs + 0.5),
     xlab = "", 
     ylab = "", 
     yaxt = "n",                 # Suppress default y-axis
     xaxt='n',
     bty = "n",                  # Remove the box around the plot
     main = "")
mtext(side=1, line=1.6, text='fold enrichment', cex=0.8)

xats = rep(1:9, 4) * rep(10^(-2:1),each=9)
xbigs = 10^(-2:1)
axis(side=1, at=xats, tck=-0.025, lwd.ticks=1, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd.ticks=1, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbigs)

# 2. Add a vertical reference line at 1
abline(v = 1, lty = 2, col = "gray60")

# 3. Plot EYH Group (#FE7F03)
segments(x0 = results$eyh_l95, y0 = y_pos + offset, 
         x1 = results$eyh_u95, y1 = y_pos + offset, 
         col = color_eyh, lwd = 2)
points(results$eyh_es, y_pos + offset, pch = 19, col = color_eyh)
mtext(side=4, line=0.25, text=ifelse(results$eyh_bonf < 0.05, "*", ""), at=y_pos + offset)

# 4. Plot LCZ Group (#7B01FC)
segments(x0 = results$lcz_l95, y0 = y_pos - offset, 
         x1 = results$lcz_u95, y1 = y_pos - offset, 
         col = color_lcz, lwd = 2)
points(results$lcz_es, y_pos - offset, pch = 19, col = color_lcz)
mtext(side=4, line=0.25, text=ifelse(results$lcz_bonf < 0.05, "*", ""), at=y_pos - offset)

# 5. Add Y-axis labels with the vertical line removed
axis(2, 
     at = y_pos, 
     labels = results$loc, 
     las = 1, 
     lwd = 0,                    # Remove the axis line
     lwd.ticks = 0,              # Remove the tick lines
     cex.axis = 0.8)

# 6. Add legend
par(xpd=T)
legend(x=3, y=0,
       legend = c("EYH", "LCZ"), 
       col = c(color_eyh, color_lcz), 
       pch = 19, 
       bty = "n")
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1


lcz_allcano %>%
  filter(p_ebm < 0.05 / nrow(lcz_allcano)) %>%
  inner_join(locs, by='gene') %>%
  filter(both %in% c("Secreted","Cell membrane")) -> lcz_secreted_membrane_hits




# panel E 


# 1. Load the datasets
eyh_full = read_tsv("output/gsea_eyh.tsv", col_types = cols())
lcz_full = read_tsv("output/gsea_lcz.tsv", col_types = cols())

# 2. Define cleaning function
clean_pathway <- function(p) {
  p %>% 
    tolower() %>% 
    str_replace("^gobp_", "") %>% 
    str_replace_all("_", " ") %>%
    str_replace_all("l amino","L amino")
}

# 3. Process and find Top 10 for each
eyh_full <- eyh_full %>% mutate(clean = clean_pathway(pathway))
lcz_full <- lcz_full %>% mutate(clean = clean_pathway(pathway))

eyh_top10 <- eyh_full %>% 
  arrange(padj) %>% 
  head(10) %>% 
  arrange(NES >= 0, padj) # Sort within block: negatives first, then by significance

lcz_top10 <- lcz_full %>% 
  arrange(padj) %>% 
  head(10) %>% 
  arrange(NES >= 0, padj) # Sort within block: negatives first, then by significance

# 4. Join datasets and define row order
# We put LCZ hits at the bottom and EYH hits at the top
shared_paths = intersect(lcz_top10$clean, eyh_top10$clean)
lcz_only_paths = setdiff(lcz_top10$clean, eyh_top10$clean)
eyh_only_paths = setdiff(eyh_top10$clean, lcz_top10$clean)

# Final sequence from bottom of plot to top of plot
ordered_paths = c(lcz_only_paths, shared_paths, eyh_only_paths)

df_plot = data.frame(clean = ordered_paths) %>%
  left_join(eyh_full %>% select(clean, nes_eyh = NES, padj_eyh = padj), by = "clean") %>%
  left_join(lcz_full %>% select(clean, nes_lcz = NES, padj_lcz = padj), by = "clean")

# Identify indices for brackets
lcz_idx_range = c(1, length(lcz_only_paths) + length(shared_paths))
eyh_idx_range = c(length(lcz_only_paths) + 1, nrow(df_plot))

# 5. Base R Plotting
par(mar = c(4, 14, 3, 6)) # Adjust margins for labels
xlims <- c(-3, 3.5)

plot(NA, NA, 
     xlim = xlims, ylim = c(0.5, nrow(df_plot) + 0.5),
     yaxt = "n", xlab = "", ylab = "",
     main = "", bty = "n", axes = FALSE)
mtext(side=1, line=2.5, text='Normalized Enrichment Score (NES)', cex=0.8)

# Custom X-axis
axis(side = 1, at = -3:3)

# Horizontal guide lines
abline(h = 1:nrow(df_plot), col = "gray95", lty = 1)

# Pathway labels
axis(2, at = 1:nrow(df_plot), labels = df_plot$clean, las = 1, cex.axis = 0.7, lwd=0, line=-0.5)

# 6. Add Dots (using the columns from the joined dataframe)
y_idx <- 1:nrow(df_plot)

# LCZ Dots (Purple)
points(df_plot$nes_lcz, y_idx, 
       col = "#7B01FC", cex = 1.6, lwd = 1.5,
       pch = ifelse(df_plot$padj_lcz < 0.05, 19, 1))

# EYH Dots (Orange)
points(df_plot$nes_eyh, y_idx, 
       col = "#FE7F03", cex = 1.6, lwd = 1.5,
       pch = ifelse(df_plot$padj_eyh < 0.05, 19, 1))

# 7. brackets
# Identify indices for brackets
lcz_idx_range = c(1, length(lcz_only_paths) + length(shared_paths))
eyh_idx_range = c(length(lcz_only_paths) + 1, nrow(df_plot))
draw_bracket = function(y_range, label, x_pos = 3.7) {
  # Vertical line
  segments(x0 = x_pos, y0 = y_range[1], x1 = x_pos, y1 = y_range[2], lwd = 1.5)
  # Top and bottom ticks
  segments(x0 = x_pos - 0.05, y0 = y_range[1], x1 = x_pos, y1 = y_range[1], lwd = 1.5)
  segments(x0 = x_pos - 0.05, y0 = y_range[2], x1 = x_pos, y1 = y_range[2], lwd = 1.5)
  # Label
  text(x = x_pos + 0.6, y = mean(y_range), labels = label, srt = 270, cex = 0.8)
}

par(xpd=T)
draw_bracket(lcz_idx_range, "Top 10 LCZ hits")
draw_bracket(eyh_idx_range, "Top 10 EYH hits")
par(xpd=F)


# 8. Legend positioning
par(xpd = TRUE)
legend(x = 3.5, y = 0, legend = c("EYH", "LCZ"), 
       col = c("#FE7F03", "#7B01FC"), pch = 19, bty = "n")

# Optional: Significance legend
legend(x = 2.5, y = -4, legend = c("padj < 0.05", "padj \u2265 0.05"), 
       pch = c(19, 1), bty = "n", horiz = TRUE, cex = 0.7)
par(xpd = FALSE)

mtext(LETTERS[panel], side=3, cex=1.5, adj = -0.1, line = 0.5); panel = panel + 1

dev.off()


