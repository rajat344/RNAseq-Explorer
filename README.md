# RNA-seq Explorer

An RNA-seq downstream analysis project focused on differential gene expression, functional enrichment, and pathway discovery across treatment conditions.

## Objective

The project analyzes publicly available human bulk RNA-seq data to identify treatment-associated changes in gene expression and investigate the biological processes associated with those changes.

## Dataset

GEO accession: GSE335921

The dataset contains 36 human stromal RNA-seq samples across three experimental microenvironment conditions and multiple treatment groups.

## Workflow

1. Dataset preparation and validation
2. Sample and experimental-design validation
3. Expression quality assessment
4. Differential expression analysis
5. Functional enrichment
6. Pathway analysis
7. Interactive visualization


### Quality Control

Initial quality-control analysis was performed using Python and pandas.

The dataset was assessed for:

- library size variation
- genes detected per sample
- low-count gene prevalence
- sample-level PCA
- sample-to-sample expression correlation

The dataset contains 36 samples and 58,735 genes.

PCA showed that the experimental microenvironment accounts for a substantial proportion of expression variation, with PC1 explaining 50.42% and PC2 explaining 23.28% of the variance.

Sample-level correlation analysis showed high overall similarity across the dataset, with a mean pairwise Pearson correlation of 0.949.

QC outputs are stored in `results/qc/` and visualizations are stored in `figures/qc/`.

## Repository structure

```text
data/          Processed expression data
metadata/      Sample and experimental metadata
scripts/       Analysis scripts
results/       Analysis outputs
figures/       Visualizations
app/           Interactive application
docs/          Project documentation