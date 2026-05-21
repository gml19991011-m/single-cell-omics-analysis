# ============================================================
# 06_visualization.R
# Generate common visualization outputs for single-cell analysis
# ============================================================

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ------------------------------------------------------------
# 1. Load annotated Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/05_annotated_seurat_object.rds")

# ------------------------------------------------------------
# 2. UMAP by cluster
# ------------------------------------------------------------

umap_by_cluster <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("UMAP by Seurat Clusters")

ggsave(
  filename = "results/figures/final_umap_by_cluster.pdf",
  plot = umap_by_cluster,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 3. UMAP by annotated cell type
# ------------------------------------------------------------

umap_by_celltype <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE
) +
  ggtitle("UMAP by Annotated Cell Type")

ggsave(
  filename = "results/figures/final_umap_by_cell_type.pdf",
  plot = umap_by_celltype,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 4. FeaturePlot for selected marker genes
# ------------------------------------------------------------

selected_genes <- c(
  "RBFOX3",
  "ELAVL4",
  "PHOX2B",
  "RET",
  "SOX10",
  "PLP1",
  "CHAT",
  "NOS1",
  "VIP",
  "TH"
)

# Keep only genes present in the dataset
selected_genes_present <- selected_genes[
  selected_genes %in% rownames(seurat_obj)
]

if (length(selected_genes_present) > 0) {
  
  feature_plot <- FeaturePlot(
    object = seurat_obj,
    features = selected_genes_present,
    reduction = "umap",
    ncol = 3
  )
  
  ggsave(
    filename = "results/figures/featureplot_selected_genes.pdf",
    plot = feature_plot,
    width = 12,
    height = 10
  )
  
} else {
  cat("No selected marker genes were found in the dataset.\n")
}

# ------------------------------------------------------------
# 5. ViolinPlot for selected marker genes
# ------------------------------------------------------------

if (length(selected_genes_present) > 0) {
  
  violin_plot <- VlnPlot(
    object = seurat_obj,
    features = selected_genes_present,
    group.by = "cell_type",
    pt.size = 0
  ) +
    NoLegend()
  
  ggsave(
    filename = "results/figures/vlnplot_selected_genes.pdf",
    plot = violin_plot,
    width = 14,
    height = 10
  )
}

# ------------------------------------------------------------
# 6. DotPlot for selected marker genes
# ------------------------------------------------------------

if (length(selected_genes_present) > 0) {
  
  dot_plot <- DotPlot(
    object = seurat_obj,
    features = selected_genes_present,
    group.by = "cell_type"
  ) +
    RotatedAxis() +
    ggtitle("Selected Marker Gene Expression")
  
  ggsave(
    filename = "results/figures/dotplot_selected_genes.pdf",
    plot = dot_plot,
    width = 10,
    height = 6
  )
}

# ------------------------------------------------------------
# 7. Cell number and proportion table
# ------------------------------------------------------------

cell_count_table <- seurat_obj@meta.data %>%
  count(cell_type, name = "cell_number") %>%
  mutate(
    proportion = cell_number / sum(cell_number)
  )

write.csv(
  cell_count_table,
  file = "results/tables/cell_type_proportions.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Bar plot of cell type proportions
# ------------------------------------------------------------

celltype_barplot <- ggplot(
  cell_count_table,
  aes(
    x = cell_type,
    y = proportion
  )
) +
  geom_col() +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Cell Type Proportions",
    x = "Cell Type",
    y = "Proportion"
  )

ggsave(
  filename = "results/figures/cell_type_proportions.pdf",
  plot = celltype_barplot,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 9. Save final metadata
# ------------------------------------------------------------

final_metadata <- seurat_obj@meta.data

write.csv(
  final_metadata,
  file = "results/tables/final_cell_metadata.csv"
)

# ------------------------------------------------------------
# 10. Print summary
# ------------------------------------------------------------

cat("Visualization completed.\n")
cat("Figures saved to results/figures/\n")
cat("Tables saved to results/tables/\n")