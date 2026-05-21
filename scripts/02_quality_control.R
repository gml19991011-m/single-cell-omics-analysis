# ============================================================
# 02_quality_control.R
# Quality control for single-cell RNA-seq data
# ============================================================

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Load raw Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/01_raw_seurat_object.rds")

# ------------------------------------------------------------
# 2. Visualize QC metrics before filtering
# ------------------------------------------------------------

qc_violin_before <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave(
  filename = "results/figures/qc_violin_before_filtering.pdf",
  plot = qc_violin_before,
  width = 12,
  height = 5
)

# ------------------------------------------------------------
# 3. Scatter plots for QC relationships
# ------------------------------------------------------------

scatter_mt <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

scatter_features <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  filename = "results/figures/qc_scatter_before_filtering.pdf",
  plot = scatter_mt + scatter_features,
  width = 10,
  height = 5
)

# ------------------------------------------------------------
# 4. Filter low-quality cells
# ------------------------------------------------------------

# These thresholds are commonly used starting points.
# They should be adjusted according to tissue type, species,
# sequencing depth, and sample quality.

seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

# ------------------------------------------------------------
# 5. Visualize QC metrics after filtering
# ------------------------------------------------------------

qc_violin_after <- VlnPlot(
  seurat_obj_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

ggsave(
  filename = "results/figures/qc_violin_after_filtering.pdf",
  plot = qc_violin_after,
  width = 12,
  height = 5
)

# ------------------------------------------------------------
# 6. Print summary information
# ------------------------------------------------------------

cat("Before filtering:\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n\n")

cat("After filtering:\n")
cat("Number of cells:", ncol(seurat_obj_filtered), "\n")
cat("Number of genes:", nrow(seurat_obj_filtered), "\n")

# ------------------------------------------------------------
# 7. Save filtered Seurat object
# ------------------------------------------------------------

saveRDS(
  seurat_obj_filtered,
  file = "data/processed/02_filtered_seurat_object.rds"
)

cat("Filtered Seurat object saved to data/processed/02_filtered_seurat_object.rds\n")