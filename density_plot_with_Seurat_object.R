# Load required libraries
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)

# Extract UMAP embeddings from Seurat object
embedding_coordinates <- seurat_object@reductions$umap@cell.embeddings %>% as.data.frame()

# Define the cells of interest based on condition
cells_of_interest <- Cells(seurat_object)[seurat_object$condition == "condition_of_interest"]

# Create the base scatter plot with light gray points
density_plot <- ggscatter(
  embedding_coordinates, 
  x = "UMAP_1", 
  y = "UMAP_2", 
  color = "lightgray", 
  size = 1
)

# Add the density plot overlay for the condition of interest
density_plot <- density_plot + 
  stat_density_2d(
    data = embedding_coordinates[cells_of_interest, ], 
    aes(fill = ..level..), 
    geom = "polygon", 
    bins = 50, 
    alpha = 0.1
  ) +
  scale_fill_gradientn(colors = c("white", '#d1e5f0', '#92c5de', '#4393c3')) +
  theme(legend.position = "right")

# Display the plot
print(density_plot)
