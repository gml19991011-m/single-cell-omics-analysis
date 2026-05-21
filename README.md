# Single-cell Omics Analysis

This repository contains a basic single-cell RNA-seq analysis workflow using R and Seurat.

The workflow includes:

1. Loading single-cell RNA-seq data
2. Quality control
3. Normalization and feature selection
4. Dimensionality reduction
5. Clustering
6. Marker gene identification
7. Manual cell type annotation
8. Visualization
9. Gene marker analysis

## Repository Structure

```text
single-cell-omics-analysis/
├── data/
│   ├── raw/              # Raw input data
│   └── processed/        # Processed Seurat objects
├── scripts/              # R scripts for single-cell analysis
├── results/
│   ├── figures/          # Output figures
│   └── tables/           # Output tables
└── docs/                 # Workflow notes
```

## Analysis Workflow

Run the scripts in the following order:

```r
source("scripts/01_load_data.R")
source("scripts/02_quality_control.R")
source("scripts/03_normalization_clustering.R")
source("scripts/04_marker_genes.R")
source("scripts/05_cell_annotation.R")
source("scripts/06_visualization.R")
source("scripts/07_marker_gene_visualization.R")
```

## Main Tools

- R
- Seurat
- dplyr
- ggplot2
- patchwork

## Research Context

This repository was built as a self-directed training project in single-cell transcriptomics analysis.

The workflow is particularly relevant to studies of neural cell diversity, enteric nervous system biology, and molecularly defined neuronal subtypes.

## Notes

This repository is intended for educational and research demonstration purposes.

The workflow should be adapted according to the specific dataset, tissue type, sequencing platform, and research question.