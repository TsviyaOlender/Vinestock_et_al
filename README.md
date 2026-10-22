# A Pre-established Fibrocartilage Population Differentiates through Successive Transcriptional Programs during Enthesis Mineralization

This repository contains the R scripts used for the single-cell RNA sequencing (scRNA-seq) analyses presented in the manuscript:

**Vinestock, R. C., Felsenthal, N., Olender, T., Krief, S., Rubin, S., Sui, X.-M., Kurtaliaj, I., Golman, M., Grinstein, M., Shahid, S., Rotkopf, R., Galloway, J. L., Foster, B. L., Thomopoulos, S., and Zelzer, E. (2026).**

*A Pre-established Fibrocartilage Population Differentiates through Successive Transcriptional Programs during Enthesis Mineralization.*

**iScience** (2026; forthcoming).

## Overview

This study investigates the cellular and transcriptional programs underlying fibrocartilage differentiation and mineralization at the tendon–bone enthesis.

Using single-cell RNA sequencing of the developing mouse Achilles enthesis, we characterized distinct cell populations and identified tenogenic, chondrogenic, and mineralization-associated transcriptional programs. The analyses provide insight into the successive transcriptional states associated with the transition from unmineralized fibrocartilage (UFC) to mineralized fibrocartilage (MFC).

The R scripts provided here support the computational analyses described in the manuscript, including cell clustering and annotation, gene module scoring, gene–gene correlation analysis, cell–cell communication inference, and comparison with a previously published scRNA-seq dataset.

## Repository Contents

### 1. `seurat_analysis_v1.R`

Main script for processing and analyzing scRNA-seq data using the Seurat R package.

The analysis includes quality-control filtering, normalization, integration of datasets, dimensionality reduction (PCA and UMAP), clustering, identification of differentially expressed genes, and automated cell-type annotation.

### 2. `score_tendon_cartilage.R`

Calculates tendon- and cartilage-associated gene module scores using predefined gene lists based on the literature.

The scores are used to characterize the transcriptional identities of cell populations and visualize tendon- and cartilage-associated expression patterns.

**Related figure:** Figure 2C.

### 3. `modules_pearson_corr.R`

Calculates pairwise Pearson correlation coefficients between selected highly variable genes using their expression profiles in the scRNA-seq dataset.

The resulting gene–gene correlation matrix is visualized as a heatmap to identify transcriptional modules associated with tenogenic, chondrogenic, and mineralization-related programs.

**Related figure:** Figure 3B.

### 4. `run_scores_fig_3D.R`

Calculates gene module scores for tenogenic, chondrogenic, and mineralization-associated programs using predefined gene sets.

The script summarizes the scores across tendon-enthesis (Tnd.Ent), unmineralized fibrocartilage (UFC), and mineralized fibrocartilage (MFC) cell populations and generates plots showing the mean scores and standard deviations.

**Related figure:** Figure 3D.

### 5. `run_cellchat_analysis.R`

Performs cell–cell communication analysis using the CellChat R package.

The analysis infers ligand–receptor interactions between cell populations, with particular emphasis on WNT signaling, and generates visualizations of predicted intercellular communication networks and signaling roles.

**Related figure:** Supplementary Figure 11A.

### 6. `analyze_cell_stem.R`

Analyzes the previously published scRNA-seq dataset from Fang et al. (2022).

The script includes data preprocessing, integration, clustering, cell-type annotation, and Seurat-based label transfer to compare cell populations between the published dataset and the Achilles enthesis dataset generated in this study.

**Related figures:** Supplementary Figures 7–8.

### 7. `function_enrichR.R`

Contains functions for automated cell-type annotation and enrichment analysis using the enrichR R package.

The functions perform gene ontology (GO), pathway, and cell-type enrichment analyses. Cell-type annotation is based primarily on the PanglaoDB database.

This script is used as a supporting function by the main Seurat analysis scripts.

## Software and Dependencies

The analyses were performed in R using packages including:

- **Seurat** – scRNA-seq processing, integration, clustering, visualization, and module scoring.
- **CellChat** – inference and visualization of cell–cell communication.
- **enrichR** – gene set enrichment and automated cell-type annotation.
- **ggplot2**, **pheatmap**, **patchwork**, and related packages – data visualization.
- **dplyr** – data manipulation.

Additional package dependencies are specified in the individual scripts.

## Usage Notes

The scripts are provided to document the computational analyses performed in the study and facilitate reuse of the analysis methods.

Some scripts require input files, including scRNA-seq count matrices, processed Seurat objects (`.rds`), and predefined gene lists. File paths in the scripts reflect the original computational environment and should be updated before running the analyses on another system.

The scripts are not organized as a fully automated end-to-end pipeline.

## Citation

If you use the code provided in this repository, please cite:

Vinestock, R. C., Felsenthal, N., Olender, T., et al. (2026). *A Pre-established Fibrocartilage Population Differentiates through Successive Transcriptional Programs during Enthesis Mineralization.* iScience.

The full publication reference and DOI will be added when available.
