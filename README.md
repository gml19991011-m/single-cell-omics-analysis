# Single-cell Omics Analysis

This repository contains a modular single-cell RNA-seq analysis workflow using **R** and **Seurat**.

The workflow includes data loading, quality control, normalization, dimensionality reduction, clustering, marker gene identification, manual cell type annotation, visualization, publication-style figure generation, and final summary table export.

## Repository Structure

```text
single-cell-omics-analysis/
├── README.md
├── LICENSE
├── .gitignore
├── data/
│   ├── raw/                  # Raw 10x Genomics input data
│   └── processed/            # Processed Seurat objects
├── scripts/                  # R scripts for analysis
├── results/
│   ├── figures/              # Standard analysis figures
│   ├── tables/               # Standard analysis tables
│   └── reports/              # Text summary reports
└── docs/                     # Workflow notes
```

## Scripts

The workflow is organized into 11 scripts.

```text
scripts/
├── 00_run_pre_annotation.R
├── 00_run_post_annotation.R
├── 01_load_data.R
├── 02_quality_control.R
├── 03_normalization_clustering.R
├── 04_marker_genes.R
├── 05_cell_annotation.R
├── 06_visualization.R
├── 07_marker_gene_visualization.R
├── 08_publication_figures.R
└── 09_annotation_and_summary.R
```

## Workflow Overview

The pipeline is divided into three stages.

---

## Stage 1. Pre-annotation analysis

Run the automatic analysis steps before manual annotation:

```r
source("scripts/00_run_pre_annotation.R")
```

This script runs:

```r
source("scripts/01_load_data.R")
source("scripts/02_quality_control.R")
source("scripts/03_normalization_clustering.R")
source("scripts/04_marker_genes.R")
```

Outputs from this stage include:

```text
results/figures/umap_clusters.pdf
results/figures/top_marker_dotplot.pdf
results/figures/top_marker_heatmap.pdf
results/tables/top10_cluster_markers.csv
results/tables/all_cluster_markers.csv
```

These files should be inspected before assigning cell type identities.

---

## Stage 2. Manual cell type annotation

After inspecting marker genes, manually update:

```text
scripts/05_cell_annotation.R
```

Specifically, modify the cluster-to-cell-type annotation section. For example:

```r
new_cluster_ids <- c(
  "0" = "CD14+ monocytes",
  "1" = "Naive CD4+ T cells",
  "2" = "CD4+ T cells",
  "3" = "CD8+ T cells",
  "4" = "B cells",
  "5" = "Cytotoxic T cells",
  "6" = "FCGR3A+ monocytes",
  "7" = "NK cells",
  "8" = "Plasmacytoid dendritic cells",
  "9" = "Innate-like T cells"
)
```

Then run:

```r
source("scripts/05_cell_annotation.R")
```

This step saves the annotated Seurat object:

```text
data/processed/05_annotated_seurat_object.rds
```

Manual annotation should always be verified using marker gene expression and biological knowledge.

---

## Stage 3. Post-annotation visualization and summary

After manual annotation has been completed, run:

```r
source("scripts/00_run_post_annotation.R")
```

This script runs:

```r
source("scripts/06_visualization.R")
source("scripts/07_marker_gene_visualization.R")
source("scripts/08_publication_figures.R")
source("scripts/09_annotation_and_summary.R")
```

Outputs from this stage include:

```text
results/figures/publication/Figure1_publication_style.pdf
results/figures/publication/pub_umap_by_cell_type.pdf
results/figures/publication/pub_marker_dotplot.pdf
results/figures/publication/pub_cell_type_proportions.pdf
results/tables/publication/final_cell_type_summary.csv
results/tables/publication/compact_marker_summary.csv
results/reports/single_cell_analysis_summary.txt
```

## Input Data Format

This workflow is designed for 10x Genomics matrix format.

Place the input files in:

```text
data/raw/filtered_feature_bc_matrix/
```

The folder should contain:

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

If the dataset contains multiple feature types, such as **Gene Expression** and **Peaks**, the workflow keeps only the **Gene Expression** matrix for scRNA-seq analysis.

## Required R Packages

The main packages used in this workflow are:

```r
Seurat
dplyr
ggplot2
patchwork
```

Install them using:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
```

## Recommended Usage

For a new dataset, the recommended workflow is:

```r
source("scripts/00_run_pre_annotation.R")
```

Then inspect marker genes and manually update:

```text
scripts/05_cell_annotation.R
```

Then run:

```r
source("scripts/05_cell_annotation.R")
source("scripts/00_run_post_annotation.R")
```

## Main Analysis Steps

### 01. Load data

```text
scripts/01_load_data.R
```

Reads 10x Genomics data and creates a Seurat object.

### 02. Quality control

```text
scripts/02_quality_control.R
```

Calculates and visualizes QC metrics, including:

- Number of detected genes
- Total UMI counts
- Mitochondrial gene percentage

### 03. Normalization and clustering

```text
scripts/03_normalization_clustering.R
```

Performs:

- Log normalization
- Highly variable gene selection
- Scaling
- PCA
- Neighbor graph construction
- Clustering
- UMAP visualization

### 04. Marker gene identification

```text
scripts/04_marker_genes.R
```

Identifies cluster-specific marker genes and exports marker tables.

### 05. Cell type annotation

```text
scripts/05_cell_annotation.R
```

Performs manual cell type annotation based on marker gene inspection.

### 06. Visualization

```text
scripts/06_visualization.R
```

Generates standard UMAP, FeaturePlot, ViolinPlot, DotPlot, and cell proportion outputs.

### 07. Marker gene visualization

```text
scripts/07_marker_gene_visualization.R
```

Generates general marker gene visualization plots.

### 08. Publication-style figures

```text
scripts/08_publication_figures.R
```

Generates cleaner publication-style figures, including:

- UMAP by cell type
- UMAP by cluster
- Marker gene dot plot
- Cell type proportion plot
- Multi-panel figure

### 09. Annotation and summary export

```text
scripts/09_annotation_and_summary.R
```

Exports final annotation tables, marker summaries, metadata, and an analysis summary report.

## Notes

This repository is intended as a self-directed training and demonstration project for single-cell transcriptomics analysis.

The pipeline should be adapted according to:

- Species
- Tissue type
- Sequencing platform
- Data quality
- Research question
- Cell type marker references

Manual cell type annotation should not be treated as fully automatic.

## Example Outputs

Key outputs include:

```text
results/figures/umap_clusters.pdf
results/figures/top_marker_dotplot.pdf
results/figures/publication/Figure1_publication_style.pdf
results/tables/top10_cluster_markers.csv
results/tables/publication/final_cell_type_summary.csv
results/reports/single_cell_analysis_summary.txt
```

## Author

Manlin Guo

## License

This repository is licensed under the MIT License.
