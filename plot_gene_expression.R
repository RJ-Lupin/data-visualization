# Draw a UMAP plot showing the normalized expression of a specific gene of interest
library(Seurat)
library(ggplot2)

# Seurat_v4 ---------------------------------------------------------------
plot_gene_expression <- function(seurat_obj, gene, order = FALSE, raster = FALSE, pt_size = 0.001, legend.position = "none") {
  if (!gene %in% rownames(seurat_obj)) stop(paste("Gene", gene, "not found."))
  max_value <- max(seurat_obj[["RNA"]]@data[gene, ])
  FeaturePlot(seurat_obj, features = gene, order = order, raster = raster, reduction = "umap", pt.size = pt_size) +
    scale_colour_gradientn(colours = rev(c("#300000", "red", "#eeeeee")),
                           breaks = c(0, max_value),
                           labels = c(0, ceiling(max_value))) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = legend.position) +
    ggtitle(gene)
}


# Seurat_v5 ---------------------------------------------------------------
plot_gene_expression <- function(seurat_obj, gene, order = FALSE, raster = FALSE, pt_size = 0.001, legend.position = "none") {
  if (!gene %in% rownames(seurat_obj)) stop(paste("Gene", gene, "not found."))
  max_value <- max(LayerData(seurat_obj.qc, assay = "RNA", layer = "data")[gene,]) # This line is different from Seurat v4
  FeaturePlot(seurat_obj, features = gene, order = order, raster = raster, reduction = "umap", pt.size = pt_size) +
    scale_colour_gradientn(colours = rev(c("#300000", "red", "#eeeeee")),
                           breaks = c(0, max_value),
                           labels = c(0, ceiling(max_value))) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = legend.position) +
    ggtitle(gene)
}
