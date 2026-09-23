# RNA-seq Explorer

A reproducible downstream RNA-seq analysis workflow for identifying treatment-associated changes in gene expression and biological pathways across different tumor microenvironment conditions.

The project combines Python-based quality control with R/DESeq2 differential expression analysis, Gene Ontology and Reactome enrichment, publication-style visualizations, and an interactive Shiny application.

---

## Overview

RNA-seq experiments can produce thousands of genes whose expression changes across experimental conditions. The challenge is to move from a raw expression matrix to statistically supported gene-level changes and biologically interpretable pathways.

This project implements a complete **downstream analysis workflow** for a public human bulk RNA-seq dataset:

```text
Processed count matrix
        |
        v
Data validation
        |
        v
Quality control
        |
        +--> Library size
        +--> Genes detected
        +--> Low-count assessment
        +--> PCA
        +--> Sample correlation
        |
        v
Differential expression
        |
        +--> DESeq2
        +--> Donor-aware experimental design
        +--> Treatment vs Control
        |
        v
Functional enrichment
        |
        +--> Gene Ontology
        +--> Biological Process
        |
        v
Pathway analysis
        |
        +--> Reactome
        |
        v
Visualization
        |
        v
Interactive Shiny Explorer
