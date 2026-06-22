ci_alpha <- 0.7

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
    dplyr::summarize(.groups='keep',
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
    mutate(x = case_when(crosshairs ~ x_mean, 
                         TRUE ~ x)) %>%
    arrange(x) -> tbl_smry
  
  if (crosshairs) {
    tbl_smry$x = tbl_smry$x_mean
  }
  
  par(mar=mar)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log=log)
  
  
  if (crosshairs) {
    segments(x0=tbl_smry$x_l95, x1=tbl_smry$x_u95, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    segments(x0=tbl_smry$x_mean, y0=tbl_smry$l95, y1=tbl_smry$u95, col=tbl_smry$color, lwd=1.5)
  }
  
  if (!is.na(boxwidth)) {
    rect(xleft=tbl_smry$x-boxwidth, xright=tbl_smry$x+boxwidth, ybottom=tbl_smry$q25, ytop=tbl_smry$q75, border=tbl_smry$color, lwd=1.5, col=NA)
    segments(x0=tbl_smry$x-boxwidth, x1=tbl_smry$x+boxwidth, y0=tbl_smry$median, col=tbl_smry$color, lwd=1.5)
  }
  
  set.seed(randseed)
  if (!is.null(pchout)) {
    points(x=jitter(tbl$x,amount=jitamt), y=tbl$y, col=pchout, pch=tbl$pch, bg=pchbg)
  } else {
    #points(x=jitter(tbl$x,amount=jitamt), y=tbl$y, col=alpha(tbl$color,ci_alpha), pch=tbl$pch, bg=pchbg) 
    points(x=jitter(tbl$x,amount=jitamt), y=tbl$y, col=tbl$color, pch=tbl$pch, bg=pchbg) 
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
  
  if (!is.na(barwidth)) {
    if (bartype=='segment') {
      segments(x0=tbl_smry$x-barwidth, x1=tbl_smry$x+barwidth, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    } else if (bartype=='bar') { #add a border
      rect(xleft=tbl_smry$x-barwidth, xright=tbl_smry$x+barwidth, ybottom=rep(0,nrow(tbl_smry)), ytop=tbl_smry$mean, col= alpha(tbl_smry$color,ci_alpha), border=tbl_smry$color)
    } #change color of arrows to black
    #arrows(x0=tbl_smry$x, y0=tbl_smry$l95, y1=tbl_smry$u95, code=3, angle=90, length=arrowlength, col=alpha("black",0.8), lwd=1.5)
    #change color to same color as the datapoints
    arrows(x0=tbl_smry$x, y0=tbl_smry$l95, y1=tbl_smry$u95, code=3, angle=90, length=arrowlength, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(test)) {
    testfun = get(test) # e.g. ks.test
    control_color = control_group
    
    tbl_smry$p = as.numeric(NA)
    for (i in 1:nrow(tbl_smry)) {
      this_x = tbl_smry$x[i]
      this_rows = tbl$x == this_x & tbl$color == tbl_smry$color[i]
      ctrl_rows =  tbl$color == control_color
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
  
  return(tbl_smry)
  
}




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
              log='', test='ks.test', control_group=dmso_col, pchvar="my_pch",  mar=c(1.5,4,1,1))

abline(h = 100, lty = 2, lwd = 1.5) 
dev.off()