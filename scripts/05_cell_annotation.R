# ============================================================
# 05_cell_annotation.R
# Manual cell type annotation based on marker genes
# ============================================================

# Load required packages
library(Seurat)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load clustered Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/03_clustered_seurat_object.rds")

# ------------------------------------------------------------
# 2. Check cluster identities
# ------------------------------------------------------------

cat("Current cluster identities:\n")
print(levels(seurat_obj))

# ------------------------------------------------------------
# 3. Define canonical marker genes
# ------------------------------------------------------------

# These markers are examples and should be adjusted according to
# tissue type, species, developmental stage, and research question.

canonical_markers <- c(
  # Neural and enteric nervous system-related markers
  "RBFOX3",   # mature neurons
  "ELAVL4",   # neuronal marker
  "PHOX2B",   # enteric neuron / autonomic lineage
  "RET",      # enteric neural crest / ENS lineage
  "SOX10",    # neural crest / glial lineage
  "PLP1",     # glial cells
  
  # Enteric neuronal subtype-related markers
  "CHAT",     # cholinergic neurons
  "SLC18A3",  # vesicular acetylcholine transporter
  "NOS1",     # nitrergic neurons
  "VIP",      # vasoactive intestinal peptide
  "TH",       # catecholaminergic neurons
  "SST",      # somatostatin-related neuronal subtype
  "CALB2",    # calretinin
  
  # Non-neural cell markers
  "EPCAM",    # epithelial cells
  "COL1A1",   # fibroblasts
  "PECAM1",   # endothelial cells
  "PTPRC",    # immune cells
  "CD3D",     # T cells
  "MS4A1",    # B cells
  "LYZ"       # myeloid cells
)

# Keep only markers present in the dataset
canonical_markers_present <- canonical_markers[
  canonical_markers %in% rownames(seurat_obj)
]

cat("Number of marker genes found in dataset:",
    length(canonical_markers_present), "\n")

# ------------------------------------------------------------
# 4. Visualize canonical marker expression
# ------------------------------------------------------------

marker_dotplot <- DotPlot(
  object = seurat_obj,
  features = canonical_markers_present
) +
  RotatedAxis()

ggsave(
  filename = "results/figures/canonical_marker_dotplot.pdf",
  plot = marker_dotplot,
  width = 12,
  height = 6
)

# Feature plots for selected markers
feature_plot <- FeaturePlot(
  object = seurat_obj,
  features = canonical_markers_present,
  reduction = "umap",
  ncol = 4
)

ggsave(
  filename = "results/figures/canonical_marker_featureplots.pdf",
  plot = feature_plot,
  width = 14,
  height = 10
)

# ------------------------------------------------------------
# 5. Manual annotation template
# ------------------------------------------------------------

# IMPORTANT:
# The annotation below is only a template.
# After checking marker expression, replace Cell_type_1, Cell_type_2, etc.
# with real cell type names.

# Example:
# new_cluster_ids <- c(
#   "0" = "Enteric neurons",
#   "1" = "Enteric glia",
#   "2" = "Fibroblasts",
#   "3" = "Epithelial cells",
#   "4" = "Immune cells"
# )

current_clusters <- levels(seurat_obj)

new_cluster_ids <- setNames(
  paste0("Cluster_", current_clusters),
  current_clusters
)

# ------------------------------------------------------------
# 6. Apply annotation
# ------------------------------------------------------------

seurat_obj <- RenameIdents(
  object = seurat_obj,
  new_cluster_ids
)

seurat_obj$cell_type <- Idents(seurat_obj)

# ------------------------------------------------------------
# 7. Visualize annotated cell types
# ------------------------------------------------------------

celltype_umap <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE
)

ggsave(
  filename = "results/figures/umap_cell_annotation.pdf",
  plot = celltype_umap,
  width = 8,
  height = 6
)

# ------------------------------------------------------------
# 8. Save annotated metadata
# ------------------------------------------------------------

metadata <- seurat_obj@meta.data

write.csv(
  metadata,
  file = "results/tables/cell_metadata_with_annotation.csv"
)

# ------------------------------------------------------------
# 9. Save annotated Seurat object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file = "data/processed/05_annotated_seurat_object.rds"
)

cat("Manual cell annotation template completed.\n")
cat("Annotated Seurat object saved to data/processed/05_annotated_seurat_object.rds\n")