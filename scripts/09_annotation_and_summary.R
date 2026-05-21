# ============================================================
# 09_annotation_and_summary.R
# Export final annotation, marker summaries, and analysis report
# ============================================================

library(Seurat)
library(dplyr)

# ------------------------------------------------------------
# 1. Load annotated Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/05_annotated_seurat_object.rds")

# ------------------------------------------------------------
# 2. Define output folders
# ------------------------------------------------------------

dir.create("results/tables/publication", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reports", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Export cluster-to-cell-type annotation table
# ------------------------------------------------------------
annotation_table <- seurat_obj@meta.data %>%
  select(seurat_clusters, cell_type) %>%
  distinct() %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  arrange(as.numeric(seurat_clusters))

write.csv(
  annotation_table,
  file = "results/tables/publication/final_cluster_annotation.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Export cell type composition table
# ------------------------------------------------------------

cell_type_summary <- seurat_obj@meta.data %>%
  count(seurat_clusters, cell_type, name = "cell_number") %>%
  mutate(
    total_cells = sum(cell_number),
    proportion = cell_number / total_cells,
    percentage = proportion * 100
  ) %>%
  arrange(seurat_clusters)

write.csv(
  cell_type_summary,
  file = "results/tables/publication/final_cell_type_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Export final metadata
# ------------------------------------------------------------

final_metadata <- seurat_obj@meta.data

write.csv(
  final_metadata,
  file = "results/tables/publication/final_metadata_for_publication.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Export top marker genes per cluster
# ------------------------------------------------------------

all_markers <- read.csv("results/tables/all_cluster_markers.csv") %>%
  mutate(cluster = as.character(cluster))

top5_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

top10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top5_markers,
  file = "results/tables/publication/top5_cluster_markers_for_publication.csv",
  row.names = FALSE
)

write.csv(
  top10_markers,
  file = "results/tables/publication/top10_cluster_markers_for_publication.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Create compact marker summary table
# ------------------------------------------------------------

marker_summary <- top5_markers %>%
  group_by(cluster) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  ) %>%
  left_join(
    annotation_table,
    by = c("cluster" = "seurat_clusters")
  ) %>%
  select(cluster, cell_type, top_markers)

write.csv(
  marker_summary,
  file = "results/tables/publication/compact_marker_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Save final Seurat object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file = "data/processed/final_annotated_seurat_object.rds"
)

# ------------------------------------------------------------
# 9. Generate text report
# ------------------------------------------------------------

report_lines <- c(
  "Single-cell RNA-seq Analysis Summary",
  "====================================",
  "",
  paste0("Total number of cells: ", ncol(seurat_obj)),
  paste0("Total number of genes: ", nrow(seurat_obj)),
  paste0("Number of Seurat clusters: ", length(unique(seurat_obj$seurat_clusters))),
  paste0("Number of annotated cell types: ", length(unique(seurat_obj$cell_type))),
  "",
  "Cell type composition:",
  paste0(
    cell_type_summary$cell_type,
    " (cluster ",
    cell_type_summary$seurat_clusters,
    "): ",
    cell_type_summary$cell_number,
    " cells, ",
    round(cell_type_summary$percentage, 2),
    "%"
  ),
  "",
  "Top marker genes per cluster:",
  paste0(
    "Cluster ",
    marker_summary$cluster,
    " - ",
    marker_summary$cell_type,
    ": ",
    marker_summary$top_markers
  ),
  "",
  "This report was generated automatically by scripts/09_annotation_and_summary.R."
)

writeLines(
  report_lines,
  con = "results/reports/single_cell_analysis_summary.txt"
)

# ------------------------------------------------------------
# 10. Print summary
# ------------------------------------------------------------

cat("Final annotation and summary export completed successfully.\n")
cat("Tables saved to results/tables/publication/\n")
cat("Report saved to results/reports/single_cell_analysis_summary.txt\n")
cat("Final Seurat object saved to data/processed/final_annotated_seurat_object.rds\n")
