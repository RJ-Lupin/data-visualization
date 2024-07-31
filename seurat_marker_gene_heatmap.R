# Load required libraries
library(Seurat)
library(pheatmap)

# Define the list of genes of interest
geneset <- c("PDCD1", "HAVCR2", "TIGIT", "LAG3", "LAYN", "CTLA4", "GZMK", "GZMA", "GZMB", "PRF1", "FGFBP2", "KLRG1", "FCGR3A", "CX3CR1")

# Extract log-normalized expression data from the Seurat object
logcounts <- GetAssayData(seurat_obj, slot = "data")
logcounts <- logcounts[geneset, ]

# Scale the data
scaled_mat <- t(scale(t(logcounts)))
scaled_mat <- na.omit(scaled_mat)

# Extract cluster information
cluster_df <- data.frame(
  cluster = seurat_obj$seurat_clusters,
  row.names = colnames(seurat_obj)
)

# Function to calculate mean expression for each cluster
calculate_cluster_means <- function(mat, clusters) {
  cluster_means <- lapply(levels(clusters), function(i) {
    rowMeans(mat[, clusters == i, drop = FALSE])
  })
  do.call(rbind, cluster_means)
}

# Calculate mean expression for each cluster
cluster_means <- calculate_cluster_means(scaled_mat, cluster_df$cluster)
rownames(cluster_means) <- paste0("Cluster_", levels(cluster_df$cluster))

# Set up color palette and range for heatmap
mypalette <- colorRampPalette(c("#053061", "#2166ac", "#f7f7f7", "#d6604d", "#9e0142"))(255)
range <- 1.5

# Generate heatmap
pheatmap(
  cluster_means,
  color = mypalette,
  breaks = seq(-range, range, length.out = 255),
  na_col = "#e0e0e0",
  show_colnames = TRUE,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  cluster_rows = TRUE
)
