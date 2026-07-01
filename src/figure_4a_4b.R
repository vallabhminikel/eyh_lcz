library(flowCore)
library(ggplot2)
library(ggcyto)
library(dplyr)
library(ggridges)

###figure 4 flow diagram### 

u251mg_fcs_filenames <- list(
  "data/flow_cytometry/export_U251MG PrP neg_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP EYH 1uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP EYH 10uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP EYH 20uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP EYH 30uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP LCZ 1uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP LCZ 10uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP LCZ 20uM_Singlets.fcs",
  "data/flow_cytometry/export_U251MG PrP pos_Singlets.fcs"
)
hek_fcs_filenames <- list(
  "data/flow_cytometry/export_Tube - 1 neg ctrl_Data Source - 1_Singlets.fcs",
  "data/flow_cytometry/export_Tube - 2 pos ctrl_[5 mL Tubes] Data Source - 1_Singlets.fcs",
  "data/flow_cytometry/export_Tube - 3 50uM EYH_[5 mL Tubes] Data Source - 1_Single Cells.fcs",
  "data/flow_cytometry/export_Tube - 4 10uM EYH_[5 mL Tubes] Data Source - 1_Singlets.fcs"
)

u251mg_file_labels <- c(
  "Negative Control",
  "EYH 1 µM",
  "EYH 10 µM",
  "EYH 20 µM",
  "EYH 30 µM",
  "LCZ 1 µM",
  "LCZ 10 µM",
  "LCZ 20 µM",
  "Positive Control"
)
hek_file_labels <- c(
  "Negative Control",
  "Positive Control",
  "EYH 50 µM",
  "EYH 10 µM"
)

#extract only the specified channel values for each .fcs file
u251mg_channel <- 'FL11-A'
u251mg_fcs_values <- list()
for(i in seq_along(u251mg_fcs_filenames)) {
  u251mg_fcs <- read.FCS(u251mg_fcs_filenames[[i]], truncate_max_range = FALSE)
  u251mg_fcs_values[[i]] <- exprs(u251mg_fcs)[, u251mg_channel]
}

hek_channel <- 'FL10-A'
hek_fcs_values <- list()
for(i in seq_along(hek_fcs_filenames)) {
  hek_fcs <- read.FCS(hek_fcs_filenames[[i]], truncate_max_range = FALSE)
  hek_fcs_values[[i]] <- exprs(hek_fcs)[, hek_channel]
}

#assign the file label names to the extracted data, then generate a df for ridgeline plotting in ggplot
names(u251mg_fcs_values) <- u251mg_file_labels
df_u251mg <- bind_rows(lapply(names(u251mg_fcs_values), function(name) {
  data.frame(
    value = u251mg_fcs_values[[name]],
    sample = name
  )
}))

names(hek_fcs_values) <- hek_file_labels
df_hek <- bind_rows(lapply(names(hek_fcs_values), function(name) {
  data.frame(
    value = hek_fcs_values[[name]],
    sample = name
  )
}))

df_u251mg <- df_u251mg %>%
  mutate(sample = factor(sample, levels = c(
    "LCZ 20 µM",
    "LCZ 10 µM",
    "LCZ 1 µM",
    "EYH 30 µM",
    "EYH 20 µM",
    "EYH 10 µM",
    "EYH 1 µM",
    "Positive Control",
    "Negative Control"
  )))

df_u251mg <- df_u251mg %>%
  mutate(flow_colors = case_when(
    grepl("EYH", sample) ~ "#FE7F03",
    grepl("LCZ", sample) ~ "#7B01FC",
    sample == "Negative Control" ~ "lightgrey",
    sample == "Positive Control" ~ "darkgrey"
  ),
  transp = case_when(
    sample == "EYH 1 µM"  ~ 0.2,
    sample == "EYH 10 µM" ~ 0.4,
    sample == "EYH 20 µM" ~ 0.6,
    sample == "EYH 30 µM" ~ 0.8,
    sample == "LCZ 1 µM"  ~ 0.2,
    sample == "LCZ 10 µM" ~ 0.4,
    sample == "LCZ 20 µM" ~ 0.6,
    sample == "Negative Control" ~ 0.5,
    sample == "Positive Control" ~ 0.5
  ))

df_hek <- df_hek %>%
  mutate(sample = factor(sample, levels = c(
    "EYH 50 µM",
    "EYH 10 µM",
    "Positive Control",
    "Negative Control"
  )))


df_hek <- df_hek %>%
  mutate(flow_colors = case_when(
    grepl("EYH", sample) ~ "#FE7F03",
    #grepl("LCZ", sample) ~ "#7B01FC",
    sample == "Negative Control" ~ "lightgrey",
    sample == "Positive Control" ~ "darkgrey"
  ),
  transp = case_when(
    sample == "EYH 10 µM" ~ 0.4,
    sample == "EYH 50 µM" ~ 1,
    sample == "Negative Control" ~ 0.5,
    sample == "Positive Control" ~ 0.5
  )
  )

#set custom x-axis breaks for logicle scaling
#want to copy flowjo output, where large tick at -10^3, -100, 0, 100, 10^3, but only 0 is annotated
#x-axis is too crowded to add -100 and 100, skipping that and having axis ticks instead
#i cant add 0 directly, because it doesn't actually sit flush with the exponents, so i have to use 10^0
custom_breaks <- c(-1e3, -100,1,100, 1e3, 1e4, 1e5, 1e6, 1e7)
custom_labels <- c(
  "",
  "",
  expression(10^0),
  "",
  "",
  expression(10^4),
  expression(10^5),
  expression(10^6),
  expression(10^7)
)


common_theme_hist <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.minor.x = element_line(color = "black", linewidth = 0.25), 
    axis.ticks.length = unit(0.5, "lines"),
    axis.text.y = element_text(color = "black", size=25, hjust=0),
    axis.text.x = element_text(color = "black", size=25),
    axis.title.x = element_text(color = "black", size=25),
    axis.title.y = element_text(color = "black", size=25),
    axis.line.x = element_line(color = 'black', linewidth = 0.3),
    plot.title = element_text(hjust = 0, size=28),
    legend.text = element_text(size = 22),
    legend.key.size = unit(1.5, "lines")
  )

#now plot a ridgeline histogram
#having y-axis with histogram counts from https://ggplot2.tidyverse.org/reference/aes_eval.html
u251mg_dataviz <- ggplot(df_u251mg, aes(x = value, y = sample, fill = flow_colors, alpha = transp, height = after_stat(count))) +  
  geom_density_ridges(aes(color = flow_colors), stat = 'binline', bins = 2048, scale = 0.95) +  
  scale_fill_identity() +  
  scale_color_identity() + 
  scale_alpha_identity() + 
  scale_x_flowjo_biexp(
    widthBasis = -100,      
    maxValue   = 16777216, 
    pos        = 6.2247198959,      
    neg        = 0,
    breaks = custom_breaks, 
    labels = custom_labels,
    minor_breaks = seq(-100, 100, by = 20),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  labs(x = "AF647 (FL11-A)", y = "", title = 'U251MG') +
  scale_y_discrete(position = 'right', expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(-1e3, 1e6)) +
  common_theme_hist 
 
hek293_dataviz <- ggplot(df_hek, aes(x = value, y = sample, fill = flow_colors, alpha = transp, height = after_stat(count))) +  
  geom_density_ridges(aes(color = flow_colors), stat = 'binline', bins = 2048, scale = 0.95) +  
  scale_fill_identity() +  
  scale_color_identity() + 
  scale_alpha_identity() + 
  scale_x_flowjo_biexp(
    widthBasis = -100,      
    maxValue   = 16777216, 
    pos        = 6.2247198959,      
    neg        = 0,
    breaks = custom_breaks, 
    labels = custom_labels,
    minor_breaks = seq(-100, 100, by = 20),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  labs(x = "AF647 (FL10-A)", y = "", title = 'HEK293') +
  scale_y_discrete(position = 'right', expand = expansion(mult = c(0, 0))) +
  coord_cartesian(xlim = c(-1e3, 1e6)) +
  common_theme_hist 

if (!dir.exists("display_items")) {
  dir.create("display_items")
}

png("display_items/figure_4_a_u251mg.png", width = 6000, height = 5000, units = "px", res=600)
print(u251mg_dataviz)
dev.off()

png("display_items/figure_4_b_hek293.png", width = 6000, height = 3000, units = "px", res=600)
print(hek293_dataviz)
dev.off()













