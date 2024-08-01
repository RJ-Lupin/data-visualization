# Load required libraries
library(ggplot2)
library(ggalluvial)
library(cowplot)
library(tibble)
library(reshape)
library(tidyverse)
library(scales)

dt <- table(seurat_obj$dpi, seurat_obj$seurat_clusters)
dt.proportion <- dt / rowSums(dt)

# Convert table to data frame and reshape to long format
df.proportion <- as.data.frame.table(dt.proportion)
colnames(df.proportion) <- c("dpi", "cluster", "value")

df.proportion <- df.proportion %>%
  mutate(
    cluster = as.factor(cluster),
    dpi = factor(dpi, levels = c("dpi0", "dpi1", "dpi2", "dpi5", "dpi7", "dpi10", "dpi14"))
  )

# Define color palette
color_palette <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928', '#e0e0e0')

# Create the plot
ggplot(df.proportion, aes(y = value, x = dpi, fill = cluster)) +
  geom_flow(aes(alluvium = cluster), color = "white", curve_type = "linear", width = 0.5, alpha = 0.4, linewidth = 1) +
  geom_col(width = 0.5, color = "black") +
  scale_fill_manual(name = NULL, values = alpha(color_palette, 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  cowplot::theme_minimal_hgrid() +
  labs(title = "Cell type proportion", x = "Timepoint", y = "Relative abundance") + 
  theme(panel.grid = element_blank())
