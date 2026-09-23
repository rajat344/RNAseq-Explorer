# RNA-seq Explorer Pipeline

## Overview

RNA-seq Explorer performs downstream analysis of human bulk RNA-seq expression data.

The workflow includes:

1. Data preparation and validation
2. Quality control
3. Experimental design validation
4. Differential expression analysis
5. Functional enrichment
6. Reactome pathway analysis
7. Visualization

## Analysis workflow

### 1. Data preparation

Processed expression data and sample metadata are validated before downstream analysis.

### 2. Quality control

The dataset is assessed using:

- Library size
- Genes detected per sample
- Low-count gene prevalence
- Principal component analysis
- Sample-to-sample Pearson correlation

### 3. Differential expression

Treatment groups are compared against control conditions to identify genes showing statistically significant expression changes.

### 4. Functional enrichment

Differentially expressed genes are analyzed for enriched Gene Ontology biological processes.

### 5. Pathway analysis

ReactomePA is used to investigate biological pathways associated with the identified gene sets.

### 6. Visualization

The pipeline generates plots for:

- PCA
- Sample correlation
- Differential expression
- GO enrichment
- Reactome enrichment

## Output organization

Analysis tables are stored in:

`results/`

Generated figures are stored in:

`figures/`

Analysis scripts are stored in:

`scripts/`

## Reproducibility

The analysis is divided into separate scripts so that individual stages can be inspected and rerun independently.