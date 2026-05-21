# ============================================================
# 03_normalization_clustering.R
# Normalization, dimensionality reduction, and clustering
# ============================================================

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Load filtered Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/02_filtered_seurat_object.rds")

# ------------------------------------------------------------
# 2. Normalize data
# ------------------------------------------------------------

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ------------------------------------------------------------
# 3. Identify highly variable genes
# ------------------------------------------------------------

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# Get top 10 variable genes
top10_variable_genes <- head(VariableFeatures(seurat_obj), 10)

# Plot variable genes
variable_feature_plot <- VariableFeaturePlot(seurat_obj)

variable_feature_plot_labeled <- LabelPoints(
  plot = variable_feature_plot,
  points = top10_variable_genes,
  repel = TRUE
)

ggsave(
  filename = "results/figures/variable_features.pdf",
  plot = variable_feature_plot_labeled,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 4. Scale data
# ------------------------------------------------------------

all_genes <- rownames(seurat_obj)

seurat_obj <- ScaleData(
  seurat_obj,
  features = all_genes
)

# ------------------------------------------------------------
# 5. Run PCA
# ------------------------------------------------------------

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj)
)

# PCA plot
pca_plot <- DimPlot(
  seurat_obj,
  reduction = "pca"
)

ggsave(
  filename = "results/figures/pca_plot.pdf",
  plot = pca_plot,
  width = 7,
  height = 6
)

# Elbow plot for choosing principal components
elbow_plot <- ElbowPlot(
  seurat_obj,
  ndims = 50
)

ggsave(
  filename = "results/figures/elbow_plot.pdf",
  plot = elbow_plot,
  width = 7,
  height = 5
)

# ------------------------------------------------------------
# 6. Find neighbors and clusters
# ------------------------------------------------------------

# The number of dimensions can be adjusted after checking the elbow plot.
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:20
)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.5
)

# ------------------------------------------------------------
# 7. Run UMAP
# ------------------------------------------------------------

seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:20
)

# UMAP plot by cluster
umap_cluster_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE
)

ggsave(
  filename = "results/figures/umap_clusters.pdf",
  plot = umap_cluster_plot,
  width = 7,
  height = 6
)

# ------------------------------------------------------------
# 8. Print summary information
# ------------------------------------------------------------

cat("Normalization, PCA, UMAP, and clustering completed.\n")
cat("Number of clusters:", length(levels(seurat_obj)), "\n")

# ------------------------------------------------------------
# 9. Save clustered Seurat object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file = "data/processed/03_clustered_seurat_object.rds"
)

cat("Clustered Seurat object saved to data/processed/03_clustered_seurat_object.rds\n")