# ============================================================
# 07_marker_gene_visualization.R
# General marker gene visualization for single-cell RNA-seq data
# ============================================================

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Load annotated Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/05_annotated_seurat_object.rds")

# ------------------------------------------------------------
# 2. Define marker genes of interest
# ------------------------------------------------------------

# This is a general marker list.
# Modify this list according to the dataset, tissue type, species,
# and research question.

marker_genes <- c(
  # General neural markers
  "RBFOX3",
  "ELAVL4",
  "MAP2",
  "TUBB3",
  
  # Glial and neural crest-related markers
  "SOX10",
  "PLP1",
  "GFAP",
  
  # Epithelial and stromal markers
  "EPCAM",
  "COL1A1",
  "DCN",
  
  # Endothelial markers
  "PECAM1",
  "VWF",
  
  # Immune cell markers
  "PTPRC",
  "CD3D",
  "MS4A1",
  "LYZ",
  
  # Proliferation markers
  "MKI67",
  "TOP2A"
)

# Keep only marker genes present in the dataset
marker_genes_present <- marker_genes[
  marker_genes %in% rownames(seurat_obj)
]

cat("Number of marker genes provided:", length(marker_genes), "\n")
cat("Number of marker genes found in dataset:", length(marker_genes_present), "\n")

# ------------------------------------------------------------
# 3. Stop if no marker genes are found
# ------------------------------------------------------------

if (length(marker_genes_present) == 0) {
  stop("None of the selected marker genes were found in the dataset. Please check gene symbols and species.")
}

# ------------------------------------------------------------
# 4. FeaturePlot for marker genes
# ------------------------------------------------------------

feature_plot <- FeaturePlot(
  object = seurat_obj,
  features = marker_genes_present,
  reduction = "umap",
  ncol = 4
)

ggsave(
  filename = "results/figures/general_marker_featureplot.pdf",
  plot = feature_plot,
  width = 14,
  height = 10
)

# ------------------------------------------------------------
# 5. DotPlot for marker genes
# ------------------------------------------------------------

dot_plot <- DotPlot(
  object = seurat_obj,
  features = marker_genes_present,
  group.by = "cell_type"
) +
  RotatedAxis() +
  ggtitle("Marker Gene Expression by Cell Type")

ggsave(
  filename = "results/figures/general_marker_dotplot.pdf",
  plot = dot_plot,
  width = 12,
  height = 6
)

# ------------------------------------------------------------
# 6. ViolinPlot for marker genes
# ------------------------------------------------------------

violin_plot <- VlnPlot(
  object = seurat_obj,
  features = marker_genes_present,
  group.by = "cell_type",
  pt.size = 0,
  ncol = 4
) +
  NoLegend()

ggsave(
  filename = "results/figures/general_marker_violinplot.pdf",
  plot = violin_plot,
  width = 14,
  height = 10
)

# ------------------------------------------------------------
# 7. Average expression table
# ------------------------------------------------------------

average_expression <- AverageExpression(
  object = seurat_obj,
  features = marker_genes_present,
  group.by = "cell_type",
  assays = "RNA"
)

average_expression_matrix <- average_expression$RNA

write.csv(
  average_expression_matrix,
  file = "results/tables/general_marker_average_expression.csv"
)

# ------------------------------------------------------------
# 8. Print summary
# ------------------------------------------------------------

cat("General marker gene visualization completed.\n")
cat("Figures saved to results/figures/\n")
cat("Average expression table saved to results/tables/general_marker_average_expression.csv\n")