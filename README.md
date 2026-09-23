
# RNA-seq Explorer

A reproducible downstream RNA-seq analysis workflow for identifying treatment-associated gene expression changes and biological pathways across different tumor microenvironment conditions.

The project combines Python-based quality control with R/DESeq2 differential expression, GO/Reactome enrichment analysis, publication-style visualizations, and an interactive Shiny application.

## Project Overview

This project analyzes publicly available human bulk RNA-seq data from GEO accession **GSE335921**.

The dataset contains **36 human stromal RNA-seq samples** across multiple treatment groups and three experimental microenvironment conditions.

The main objective is to determine how different treatments affect gene expression within each microenvironment and to identify biological processes and pathways associated with the observed changes.

## Workflow

```text
Processed RNA-seq counts
        |
        v
Data preparation & validation
        |
        v
Quality control
  - Library size
  - Genes detected
  - Low-count assessment
  - PCA
  - Sample correlation
        |
        v
Differential expression
  - DESeq2
  - Donor-aware design
  - Treatment vs Control
        |
        v
Functional enrichment
  - Gene Ontology
  - Biological Process
        |
        v
Pathway analysis
  - Reactome
        |
        v
Visualization
        |
        v
Interactive Shiny Explorer
