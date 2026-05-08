# BF530 Final Project — RNA-seq Expression Analysis Shiny App
**Christine Snow | Spring 2026**

## Overview
An interactive R Shiny application for exploring RNA-seq gene expression data from a study comparing post-mortem prefrontal cortex (BA9) tissue from Huntington's Disease patients and neurologically normal controls (GSE64810, Labadorf et al. 2015).

## Dataset
- **GEO Accession:** GSE64810
- **Comparison:** 20 Huntington's Disease vs 49 neurologically normal controls
- **Tissue:** Human dorsolateral prefrontal cortex (Brodmann Area 9)
  
## Preprocessing Code
The R Markdown file used to preprocess the raw data is included in this repository as `shiny-markdow.Rmd`. It contains all code for:
- Loading and parsing the raw counts matrix and series matrix
- Building the sample metadata table (coldata)
- Filtering low-expressed genes
- Running DESeq2
- VST normalization
- Saving all three CSV files used by the Shiny app
  
## App Features

### Tab 1 — Sample Info
Upload `sample_info.csv` to explore sample metadata. Includes a summary table, sortable data table, and histogram of continuous variables.

### Tab 2 — Counts Exploration
Upload `normalized_counts.csv` to explore the expression matrix. Filter genes by variance percentile and minimum non-zero samples using reactive sliders. Includes a filter summary, scatter plot, clustered heatmap, PCA, and data table.

### Tab 3 — Differential Expression
Upload `de_results.csv` to explore DESeq2 results. Includes a sortable results table and volcano plot colored by significance (FDR < 0.05).

### Tab 4 — Network Analysis
Upload `normalized_counts.csv`, enter gene names one per line, and set a correlation threshold. Computes pairwise Pearson correlations and visualizes a gene correlation network. Includes a heatmap, network plot, and statistics table (degree, closeness, betweenness).

## Data Files
All processed data files are in the `data/` folder:
- `sample_info.csv` — sample metadata including condition, age of death, PMI, and RIN
- `normalized_counts.csv` — VST normalized counts matrix
- `de_results.csv` — DESeq2 differential expression results (HD vs Control)

## Preprocessing
Raw data was downloaded from GEO (GSE64810). Genes were filtered following the original authors' approach — removing genes where more than 50% of samples in either condition had zero counts (39,376 → 26,690 genes). DESeq2 was run comparing HD vs Control with condition as the only covariate. Counts were VST normalized for visualization.

## Requirements
```r
install.packages(c("shiny", "bslib", "tidyverse", "DT", 
                   "pheatmap", "ggrepel", "igraph", "showtext"))
```

## How to Run
```r
shiny::runApp("app.R")
```
