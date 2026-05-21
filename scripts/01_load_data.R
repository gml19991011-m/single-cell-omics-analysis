# ============================================================
# 01_load_data.R
# Load single-cell RNA-seq data and create a Seurat object
# ============================================================

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Set project information
# ------------------------------------------------------------

project_name <- "single_cell_omics_analysis"

# The input folder should contain:
# - matrix.mtx.gz
# - features.tsv.gz or genes.tsv.gz
# - barcodes.tsv.gz
#
# Example path for 10x Genomics output:
# data/raw/filtered_feature_bc_matrix

data_dir <- "data/raw/filtered_feature_bc_matrix"

# ------------------------------------------------------------
# 2. Load 10x Genomics data
# ------------------------------------------------------------

counts <- Read10X(data.dir = data_dir)

# ------------------------------------------------------------
# 3. Create Seurat object
# ------------------------------------------------------------

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = project_name,
  min.cells = 3,
  min.features = 200
)

# ------------------------------------------------------------
# 4. Calculate mitochondrial gene percentage
# ------------------------------------------------------------

# For human data, mitochondrial genes usually start with "MT-"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  object = seurat_obj,
  pattern = "^MT-"
)

# For mouse data, mitochondrial genes usually start with "mt-"
# If you are analyzing mouse data, use the following line instead:
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
#   object = seurat_obj,
#   pattern = "^mt-"
# )

# ------------------------------------------------------------
# 5. Print basic information
# ------------------------------------------------------------

print(seurat_obj)

cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")

# ------------------------------------------------------------
# 6. Save Seurat object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file = "data/processed/01_raw_seurat_object.rds"
)

cat("Raw Seurat object saved to data/processed/01_raw_seurat_object.rds\n")