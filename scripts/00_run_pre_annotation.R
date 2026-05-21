# ============================================================
# 00_run_pre_annotation.R
# Run the pre-annotation single-cell RNA-seq workflow
# ============================================================

# This script runs the analysis steps before manual cell type annotation:
# 01. Load data
# 02. Quality control
# 03. Normalization, PCA, UMAP, and clustering
# 04. Marker gene identification
#
# After running this script, inspect the marker genes and manually update:
# scripts/05_cell_annotation.R

# ------------------------------------------------------------
# 1. Set project root
# ------------------------------------------------------------

setwd("~/single-cell-omics-analysis")

cat("Current working directory:\n")
cat(getwd(), "\n\n")

# ------------------------------------------------------------
# 2. Check input data folder
# ------------------------------------------------------------

data_dir <- "data/raw/filtered_feature_bc_matrix"

if (!dir.exists(data_dir)) {
  stop(
    paste0(
      "Input data folder not found: ",
      data_dir,
      "\nPlease place the 10x matrix files in this folder."
    )
  )
}

cat("Input data folder found:\n")
cat(data_dir, "\n\n")

cat("Input files:\n")
print(list.files(data_dir))
cat("\n")

# ------------------------------------------------------------
# 3. Check required packages
# ------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      "\nPlease install them before running the workflow."
    )
  )
}

cat("All required packages are available.\n\n")

# ------------------------------------------------------------
# 4. Create output folders
# ------------------------------------------------------------

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("Output folders checked.\n\n")

# ------------------------------------------------------------
# 5. Run pre-annotation workflow
# ------------------------------------------------------------

cat("Step 01: Loading data...\n")
source("scripts/01_load_data.R")
cat("Step 01 completed.\n\n")

cat("Step 02: Quality control...\n")
source("scripts/02_quality_control.R")
cat("Step 02 completed.\n\n")

cat("Step 03: Normalization, PCA, UMAP, and clustering...\n")
source("scripts/03_normalization_clustering.R")
cat("Step 03 completed.\n\n")

cat("Step 04: Marker gene identification...\n")
source("scripts/04_marker_genes.R")
cat("Step 04 completed.\n\n")

# ------------------------------------------------------------
# 6. Final message
# ------------------------------------------------------------

cat("============================================================\n")
cat("Pre-annotation workflow completed successfully.\n")
cat("============================================================\n\n")

cat("Please inspect the following outputs before manual annotation:\n")
cat("- results/figures/umap_clusters.pdf\n")
cat("- results/figures/top_marker_dotplot.pdf\n")
cat("- results/figures/top_marker_heatmap.pdf\n")
cat("- results/tables/top10_cluster_markers.csv\n")
cat("- results/tables/all_cluster_markers.csv\n\n")

cat("Next step:\n")
cat("Manually update scripts/05_cell_annotation.R based on marker genes.\n")
