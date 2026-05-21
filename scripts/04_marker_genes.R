# ============================================================
# 04_marker_genes.R
# Identify marker genes for each cluster
# ============================================================

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load clustered Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/03_clustered_seurat_object.rds")

# ------------------------------------------------------------
# 2. Identify marker genes for all clusters
# ------------------------------------------------------------

# FindAllMarkers identifies genes that are enriched in each cluster.
# only.pos = TRUE means only positive markers are reported.

cluster_markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# ------------------------------------------------------------
# 3. Save all marker genes
# ------------------------------------------------------------

write.csv(
  cluster_markers,
  file = "results/tables/all_cluster_markers.csv",
  row.names = FALSE
)

cat("All cluster markers saved to results/tables/all_cluster_markers.csv\n")

# ------------------------------------------------------------
# 4. Extract top marker genes for each cluster
# ------------------------------------------------------------

top10_markers <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  )

write.csv(
  top10_markers,
  file = "results/tables/top10_cluster_markers.csv",
  row.names = FALSE
)

cat("Top 10 cluster markers saved to results/tables/top10_cluster_markers.csv\n")

# ------------------------------------------------------------
# 5. Dot plot of top marker genes
# ------------------------------------------------------------

top_genes <- unique(top10_markers$gene)

marker_dotplot <- DotPlot(
  object = seurat_obj,
  features = top_genes
) +
  RotatedAxis()

ggsave(
  filename = "results/figures/top_marker_dotplot.pdf",
  plot = marker_dotplot,
  width = 14,
  height = 8
)

# ------------------------------------------------------------
# 6. Heatmap of top marker genes
# ------------------------------------------------------------

marker_heatmap <- DoHeatmap(
  object = seurat_obj,
  features = top_genes
) +
  NoLegend()

ggsave(
  filename = "results/figures/top_marker_heatmap.pdf",
  plot = marker_heatmap,
  width = 14,
  height = 10
)

# ------------------------------------------------------------
# 7. Print summary information
# ------------------------------------------------------------

cat("Marker gene analysis completed.\n")
cat("Number of marker genes identified:", nrow(cluster_markers), "\n")
cat("Number of clusters:", length(unique(cluster_markers$cluster)), "\n")