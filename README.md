# BF530 Final Project — RNA-seq Expression Analysis Shiny App
**Christine Snow | Spring 2026**

## Overview
An interactive R Shiny application for exploring RNA-seq gene expression data from a study comparing post-mortem prefrontal cortex (BA9) tissue from Huntington's Disease patients and neurologically normal controls (GSE64810, Labadorf et al. 2015).

## Dataset
- **GEO Accession:** GSE64810
- **Comparison:** 20 Huntington's Disease vs 49 neurologically normal controls
- **Tissue:** Human dorsolateral prefrontal cortex (Brodmann Area 9)

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
- `sample_info.csv` — sam
