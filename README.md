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

## Repository structure

```text
data/          Processed expression data
metadata/      Sample and experimental metadata
scripts/       Analysis scripts
results/       Analysis outputs
figures/       Visualizations
app/           Interactive application
docs/          Project documentation