#### To visualize gene expression or signature scores in single-cell RNA-seq data
# Load required libraries
library(Seurat)
library(ggplot2)

# Define the function to generate a feature plot with detailed legend
seurat_feature_plot <- function(seurat_obj, feature, feature_type = "gene", reduction = "umap", legend_title = "Feature Score") {
  
  # Calculate the maximum value for the specified feature
  max_value <- switch(feature_type,
                      gene = max(seurat_obj@assays$RNA@data[feature, ], na.rm = TRUE),
                      signature_score = max(seurat_obj[[feature]], na.rm = TRUE),
                      stop("Invalid feature type. Please specify either 'gene' or 'signature_score'."))
  
  # Generate the feature plot
  feature_plot <- FeaturePlot(seurat_obj, features = feature, reduction = reduction) +
    scale_colour_gradientn(colours = c("#eeeeee", "red", "#300000"),
                           breaks = c(0, max_value),
                           labels = c(0, ceiling(max_value)),
                           name = legend_title) +
    theme_void() +
    theme(legend.position = "right")
  
  # Return the plot
  return(feature_plot)
}

# Example usage
# seurat_feature_plot(seurat_object, feature_type = "gene", feature = "gene_name")
# seurat_feature_plot(seurat_object, feature_type = "signature_score", feature = "signature_name")
