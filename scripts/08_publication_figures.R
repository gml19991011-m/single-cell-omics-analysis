# ============================================================
# 08_publication_figures.R
# Generate publication-style figures for single-cell RNA-seq data
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# ------------------------------------------------------------
# 1. Load annotated Seurat object
# ------------------------------------------------------------

seurat_obj <- readRDS("data/processed/05_annotated_seurat_object.rds")

# ------------------------------------------------------------
# 2. Define output folders
# ------------------------------------------------------------

dir.create("results/figures/publication", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables/publication", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Define a publication-style theme
# ------------------------------------------------------------

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 12, face = "bold", colour = "black"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )

# ------------------------------------------------------------
# 4. UMAP by annotated cell type
# ------------------------------------------------------------

umap_celltype_pub <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.55
) +
  labs(
    title = "UMAP of annotated cell types",
    color = "Cell type"
  ) +
  theme_pub +
  theme(
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

ggsave(
  filename = "results/figures/publication/pub_umap_by_cell_type.pdf",
  plot = umap_celltype_pub,
  width = 8,
  height = 6
)

ggsave(
  filename = "results/figures/publication/pub_umap_by_cell_type.png",
  plot = umap_celltype_pub,
  width = 8,
  height = 6,
  dpi = 600
)

# ------------------------------------------------------------
# 5. UMAP by Seurat cluster
# ------------------------------------------------------------

umap_cluster_pub <- DimPlot(
  object = seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.55
) +
  labs(
    title = "UMAP by Seurat clusters",
    color = "Cluster"
  ) +
  theme_pub

ggsave(
  filename = "results/figures/publication/pub_umap_by_cluster.pdf",
  plot = umap_cluster_pub,
  width = 8,
  height = 6
)

ggsave(
  filename = "results/figures/publication/pub_umap_by_cluster.png",
  plot = umap_cluster_pub,
  width = 8,
  height = 6,
  dpi = 600
)

# ------------------------------------------------------------
# 6. Marker gene DotPlot
# ------------------------------------------------------------

marker_genes <- c(
  "VCAN", "TREM1", "FCAR", "CSF3R",
  "LEF1", "MAL", "GATA3", "CD28",
  "CD8A", "CD8B", "MS4A1", "PAX5",
  "GZMH", "GZMK", "EOMES", "FCGR3A",
  "NCAM1", "NCR1", "CLEC4C", "IL23R"
)

marker_genes_present <- marker_genes[marker_genes %in% rownames(seurat_obj)]

dotplot_pub <- DotPlot(
  object = seurat_obj,
  features = marker_genes_present,
  group.by = "cell_type"
) +
  RotatedAxis() +
  labs(
    title = "Representative marker gene expression",
    x = "Marker genes",
    y = "Cell type"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(
  filename = "results/figures/publication/pub_marker_dotplot.pdf",
  plot = dotplot_pub,
  width = 11,
  height = 6
)

ggsave(
  filename = "results/figures/publication/pub_marker_dotplot.png",
  plot = dotplot_pub,
  width = 11,
  height = 6,
  dpi = 600
)

# ------------------------------------------------------------
# 7. Cell type proportion table and plot
# ------------------------------------------------------------

cell_count_table <- seurat_obj@meta.data %>%
  count(cell_type, name = "cell_number") %>%
  mutate(
    proportion = cell_number / sum(cell_number),
    percentage = proportion * 100
  ) %>%
  arrange(desc(proportion))

write.csv(
  cell_count_table,
  file = "results/tables/publication/pub_cell_type_proportions.csv",
  row.names = FALSE
)

barplot_pub <- ggplot(
  cell_count_table,
  aes(
    x = reorder(cell_type, proportion),
    y = percentage
  )
) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(
    title = "Cell type composition",
    x = "Cell type",
    y = "Percentage of cells (%)"
  ) +
  theme_pub

ggsave(
  filename = "results/figures/publication/pub_cell_type_proportions.pdf",
  plot = barplot_pub,
  width = 7,
  height = 5
)

ggsave(
  filename = "results/figures/publication/pub_cell_type_proportions.png",
  plot = barplot_pub,
  width = 7,
  height = 5,
  dpi = 600
)

# ------------------------------------------------------------
# 8. Multi-panel figure
# ------------------------------------------------------------

figure1 <- (umap_celltype_pub | dotplot_pub) / barplot_pub +
  plot_annotation(tag_levels = "A")

ggsave(
  filename = "results/figures/publication/Figure1_publication_style.pdf",
  plot = figure1,
  width = 14,
  height = 10
)

ggsave(
  filename = "results/figures/publication/Figure1_publication_style.png",
  plot = figure1,
  width = 14,
  height = 10,
  dpi = 600
)

# ------------------------------------------------------------
# 9. Print summary
# ------------------------------------------------------------

cat("Publication-style figures generated successfully.\n")
cat("Output folder: results/figures/publication/\n")
cat("Summary table: results/tables/publication/pub_cell_type_proportions.csv\n")
