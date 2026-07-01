library(flowCore)
library(ggplot2)
library(ggcyto)
library(dplyr)
library(patchwork)

supplemental_filenames <- list(
   "data/flow_cytometry/export_A375 neg_Single Cells.fcs",
   "data/flow_cytometry/export_A375 pos_singlets.fcs",
   "data/flow_cytometry/export_HT29 neg_singlets.fcs",
   "data/flow_cytometry/export_HT29 pos_singlets.fcs",
   "data/flow_cytometry/export_U87MG neg_singlets.fcs",
   "data/flow_cytometry/export_U87MG pos_singlets.fcs",
   "data/flow_cytometry/export_U251MG neg_Singlets I.fcs",
   "data/flow_cytometry/export_U251MG pos_Singlets I.fcs"
   )

supplemental_labels <- c(
  "A375 Negative",
  "A375 Positive",
  "HT29 Negative",
  "HT29 Positive",
  "U87MG Negative",
  "U87MG Positive",
  "U251MG Negative",
  "U251MG Positive"
)

#loop along every .fcs file, read it, and extract only the specified channel values
channel <- 'FL11-A'
fcs_values <- list()
for(i in seq_along(supplemental_filenames)) {
  fcs <- read.FCS(supplemental_filenames[[i]], truncate_max_range = FALSE)
  fcs_values[[i]] <- exprs(fcs)[, channel]
}

#assign the file label names to the extracted data, then generate a df that I further filter for cell line 
names(fcs_values) <- supplemental_labels
df_supplemental <- bind_rows(lapply(names(fcs_values), function(name) {
  data.frame(
    value = fcs_values[[name]],
    sample = name
  )
}))

#separate cell line vs status (Positive/Negative)
df_supplemental <- df_supplemental %>%
  mutate(
    cell_line = sub(" .*", "", sample),   #text before first space
    status = ifelse(grepl("Positive", sample), "Cas9-Positive", "Negative Control")
  )


#set custom x-axis breaks for biexponential scaling
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
common_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.minor.x = element_line(color = "black", linewidth = 0.25), # Makes minor ticks visible
    axis.ticks.length = unit(0.5, "lines"),
    axis.text.y = element_text(color = "black", size=25),
    axis.text.x = element_text(color = "black", size=25),
    axis.title.x = element_text(color = "black", size=25),
    axis.title.y = element_text(color = "black", size=25),
    axis.line.x = element_line(color = 'black', linewidth = 0.3),
    axis.line.y = element_line(color = 'black', linewidth = 0.3),
    plot.title = element_text(hjust = 0.5, face = "bold", size=28),
    legend.text = element_text(size = 22),
    legend.key.size = unit(1.5, "lines")
  )

custom_x_scale <- scale_x_flowjo_biexp(
  widthBasis = -100, maxValue = 16777216, pos = 6.2247198959, neg = 0,
  breaks = custom_breaks, 
  labels = custom_labels,
  minor_breaks = seq(-100, 100, by = 20),
  guide = guide_axis(minor.ticks = TRUE)
)


#plot for A375
p_a375 <- ggplot(df_supplemental %>% dplyr::filter(cell_line == "A375"), aes(x = value, fill = status)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 200) +
  scale_fill_manual(values = c("Negative Control" = "grey", "Cas9-Positive" = "cyan"), name = '') +
  custom_x_scale + 
  scale_y_continuous(expand = c(0, 0.05)) +
  coord_cartesian(xlim = c(-1.1e3, 1e6)) +
  labs(title = "A375", x = channel, y = "Count") +
  common_theme

#plot for HT29
p_ht29 <- ggplot(df_supplemental %>% dplyr::filter(cell_line == "HT29"), aes(x = value, fill = status)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 200) +
  scale_fill_manual(values = c("Negative Control" = "grey", "Cas9-Positive" = "cyan"), name = '') +
  custom_x_scale + 
  scale_y_continuous(expand = c(0, 0.05)) +
  coord_cartesian(xlim = c(-1.1e3, 1e6)) +
  labs(title = "HT29", x = channel, y = "Count") +
  common_theme

#plot for U87MG
p_u87mg <- ggplot(df_supplemental %>% dplyr::filter(cell_line == "U87MG"), aes(x = value, fill = status)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 200) +
  scale_fill_manual(values = c("Negative Control" = "grey", "Cas9-Positive" = "cyan"), name = '') +
  custom_x_scale + 
  scale_y_continuous(expand = c(0, 0.05)) +
  coord_cartesian(xlim = c(-1.1e3, 1e6)) +
  labs(title = "U87MG", x = channel, y = "Count") +
  common_theme

#plot for U251MG
p_u251mg <- ggplot(df_supplemental %>% dplyr::filter(cell_line == "U251MG"), aes(x = value, fill = status)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 200) +
  scale_fill_manual(values = c("Negative Control" = "grey", "Cas9-Positive" = "cyan"), name = '') +
  custom_x_scale +
  scale_y_continuous(expand = c(0, 0.05)) +
  coord_cartesian(xlim = c(-1.1e3, 1e6)) +
  labs(title = "U251MG", x = channel, y = "Count") +
  common_theme

#patchwork to combine into one figure
combined_plot <- (p_a375 | p_ht29) / (p_u87mg | p_u251mg) +
  plot_layout(guides = "collect")

print(combined_plot)

if (!dir.exists("display_items")) {
  dir.create("display_items")
}
ggsave("display_items/figure_s6b.png", plot = combined_plot, dpi=600, height = 7800, width = 9000, units = 'px')


