# ============================================================
# 00_run_post_annotation.R
# Run the post-annotation single-cell RNA-seq workflow
# ============================================================

# This script runs the analysis steps after manual cell type annotation:
# 06. General visualization
# 07. Marker gene visualization
# 08. Publication-style figures
# 09. Final annotation and summary tables
#
# IMPORTANT:
# Before running this script, make sure that:
# 1. scripts/05_cell_annotation.R has been manually updated.
# 2. source("scripts/05_cell_annotation.R") has been run successfully.
# 3. data/processed/05_annotated_seurat_object.rds exists.

# ------------------------------------------------------------
# 1. Set project root
# ------------------------------------------------------------

setwd("~/single-cell-omics-analysis")

cat("Current working directory:\n")
cat(getwd(), "\n\n")

# ------------------------------------------------------------
# 2. Check annotated Seurat object
# ------------------------------------------------------------

annotated_object <- "data/processed/05_annotated_seurat_object.rds"

if (!file.exists(annotated_object)) {
  stop(
    paste0(
      "Annotated Seurat object not found: ",
      annotated_object,
      "\nPlease run scripts/05_cell_annotation.R before running this workflow."
    )
  )
}

cat("Annotated Seurat object found:\n")
cat(annotated_object, "\n\n")

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

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/publication", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/publication", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reports", recursive = TRUE, showWarnings = FALSE)

cat("Output folders checked.\n\n")

# ------------------------------------------------------------
# 5. Run post-annotation workflow
# ------------------------------------------------------------

cat("Step 06: General visualization...\n")
source("scripts/06_visualization.R")
cat("Step 06 completed.\n\n")

cat("Step 07: General marker gene visualization...\n")
source("scripts/07_marker_gene_visualization.R")
cat("Step 07 completed.\n\n")

cat("Step 08: Publication-style figures...\n")
source("scripts/08_publication_figures.R")
cat("Step 08 completed.\n\n")

cat("Step 09: Final annotation and summary tables...\n")
source("scripts/09_annotation_and_summary.R")
cat("Step 09 completed.\n\n")

# ------------------------------------------------------------
# 6. Final message
# ------------------------------------------------------------

cat("============================================================\n")
cat("Post-annotation workflow completed successfully.\n")
cat("============================================================\n\n")

cat("Key output folders:\n")
cat("- results/figures/\n")
cat("- results/figures/publication/\n")
cat("- results/tables/\n")
cat("- results/tables/publication/\n")
cat("- results/reports/\n\n")

cat("Recommended files to inspect:\n")
cat("- results/figures/publication/Figure1_publication_style.pdf\n")
cat("- results/figures/publication/pub_umap_by_cell_type.pdf\n")
cat("- results/figures/publication/pub_marker_dotplot.pdf\n")
cat("- results/tables/publication/final_cell_type_summary.csv\n")
cat("- results/tables/publication/compact_marker_summary.csv\n")
cat("- results/reports/single_cell_analysis_summary.txt\n\n")

cat("To open publication figures in RStudio, run:\n")
cat('system("open results/figures/publication")\n\n')

cat("To open the summary report in RStudio, run:\n")
cat('system("open results/reports/single_cell_analysis_summary.txt")\n')
