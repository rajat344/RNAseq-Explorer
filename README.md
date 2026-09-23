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

Objective

The objective is to determine how different treatments affect gene expression within distinct experimental microenvironment conditions and to investigate the biological processes and pathways associated with those changes.

The project was also designed to practice a reproducible computational biology workflow combining:

Statistical analysis
Python and R programming
RNA-seq quality control
Differential expression
Functional enrichment
Data visualization
Reproducible project organization
Interactive biological data exploration
Dataset

Source: NCBI Gene Expression Omnibus (GEO)

GEO accession: GSE335921

The dataset contains:

36 human stromal RNA-seq samples
58,735 genes
Three experimental microenvironment conditions
Multiple treatment groups

The analysis uses processed expression/count data rather than starting from raw FASTQ files.

Analysis Workflow
1. Data Preparation and Validation

The processed count matrix and sample metadata are validated before downstream analysis.

The workflow checks:

Required input files
Gene identifiers
Duplicate gene IDs
Sample identifiers
Sample-to-metadata matching
Numeric count values
Negative counts
Integer count requirements
Metadata completeness
Sample ordering

These checks prevent common input and experimental-design errors from propagating into downstream analysis.

2. Quality Control

Initial expression-level quality control is performed using Python.

The QC workflow evaluates:

Library size variation
Number of genes detected per sample
Gene-level count distributions
Low-count prevalence
Sample-level PCA
Sample-to-sample expression correlation

The dataset contains 36 samples and 58,735 genes.

PCA showed:

PC1: 50.42% variance explained
PC2: 23.28% variance explained

The mean pairwise Pearson correlation across samples was 0.949.

These results provide an initial assessment of sample-level structure and expression similarity.

QC outputs are stored in:

results/qc/

QC visualizations are stored in:

figures/qc/
3. Differential Expression Analysis

Differential expression analysis is performed using DESeq2.

The analysis compares treatment groups against their corresponding control within each microenvironment.

The experimental design includes donor as a covariate:

~ donor + treatment

This accounts for donor-associated variation while estimating treatment effects.

The workflow generates differential expression results containing:

Base mean expression
log2 fold change
Standard error
Wald statistic
p-value
Adjusted p-value

Significance summaries are generated using an adjusted p-value threshold of:

FDR < 0.05

The analysis script also validates the input count matrix and metadata before running DESeq2.

Differential expression results are stored in:

results/differential_expression/

The main analysis script is:

scripts/differential_expression.R
4. Functional Enrichment

Significant genes are separated into:

Upregulated genes
Downregulated genes

Gene identifiers are mapped from Ensembl to Entrez identifiers using:

org.Hs.eg.db

Gene Ontology Biological Process enrichment is performed using:

clusterProfiler

Multiple-testing correction is applied to enrichment results.

The enrichment workflow generates:

Input gene counts
Successfully mapped genes
GO Biological Process terms
Reactome pathway results
Enrichment visualizations

Results are stored in:

results/enrichment/

Visualizations are stored in:

figures/enrichment/

The main analysis script is:

scripts/enrichment.R
Gene-ID mapping

Gene-ID mapping rates vary between comparisons.

This is expected in real-world public datasets because not every identifier necessarily maps successfully to the annotation database being used.

The project therefore reports both:

input_genes
mapped_genes

rather than assuming that every input gene was successfully annotated.

Interactive RNA-seq Explorer

The project includes an R Shiny application for interactively exploring enrichment results.

The application provides:

Gene-direction filtering
Enrichment summary tables
Comparison-level gene counts
Interactive visualization of analysis outputs

The application is located at:

app/app.R

Run it from the project root with:

shiny::runApp("app")

The application is intended as an exploration layer on top of the generated analysis results rather than as a replacement for the statistical workflow.

Visualizations

The project generates visual outputs for multiple stages of the analysis.

Quality Control

Examples include:

PCA
Sample correlation
Sample-level QC summaries
Gene-level QC summaries
Differential Expression

Examples include:

Differential expression summaries
Comparison-level visualizations
Functional Enrichment

Examples include:

GO Biological Process dot plots
Reactome enrichment plots
Comparison-specific enrichment figures

All generated figures are organized under:

figures/
Technologies
Programming Languages
Python
R
Python
pandas
NumPy
matplotlib
R / Bioconductor
DESeq2
clusterProfiler
org.Hs.eg.db
ReactomePA
enrichplot
ggplot2
Shiny
Data Source
NCBI GEO
GEO accession: GSE335921
Repository Structure
RNAseq-Explorer/
│
├── app/
│   └── app.R
│
├── data/
│   ├── processed/
│   └── raw/
│
├── metadata/
│   └── samples.csv
│
├── scripts/
│   ├── prepare_data.py
│   ├── qc.py
│   ├── pca.py
│   ├── correlation.py
│   ├── design.py
│   ├── differential_expression.R
│   ├── enrichment.R
│   └── de_visualization.R
│
├── results/
│   ├── qc/
│   ├── differential_expression/
│   └── enrichment/
│
├── figures/
│   ├── qc/
│   ├── differential_expression/
│   └── enrichment/
│
├── docs/
│   └── pipeline.md
│
├── .gitignore
└── README.md
Reproducibility

The project separates data, metadata, analysis scripts, results, figures, documentation, and application code.

This structure allows each stage of the workflow to be inspected independently.

The main analysis stages are implemented as separate scripts rather than a single monolithic program.

Important generated outputs are retained in the repository so that the analysis can be inspected without rerunning every step.

For a detailed description of the workflow, see:

docs/pipeline.md
Running the Project
Python analysis

From the project root:

python scripts/qc.py
python scripts/pca.py
python scripts/correlation.py

Additional preparation and visualization scripts are available under:

scripts/
R analysis

Run the differential expression workflow from the project root:

source("scripts/differential_expression.R")

Then run enrichment:

source("scripts/enrichment.R")

Visualization scripts are available under:

scripts/de_visualization.R
Shiny application

From the project root:

shiny::runApp("app")
Interpretation

The analysis is designed to identify statistical and functional patterns associated with treatment conditions.

The results should be interpreted at the level of:

Differentially expressed genes
Direction of expression change
Enriched biological processes
Enriched pathways
Differences between treatment and microenvironment conditions

Enrichment results provide biological context for the observed gene-level changes but do not by themselves establish causal mechanisms.

Limitations

This project is intentionally focused on downstream RNA-seq analysis.

The current workflow does not start from raw sequencing reads and therefore does not include:

FASTQ quality control
Adapter trimming
Read alignment
Transcript quantification
featureCounts
Salmon
STAR
HISAT2
HPC-based raw-read processing
Nextflow workflow orchestration

Therefore, this repository should not be described as a complete raw-read RNA-seq processing pipeline.

It is a downstream analysis workflow beginning with processed expression/count data.

Additional limitations include:

Gene-ID mapping rates vary between comparisons.
Functional enrichment depends on successful identifier mapping.
Results are based on a single public dataset.
Computational enrichment does not establish experimental causality.
Biological interpretation requires appropriate experimental validation.
What This Project Demonstrates

This project demonstrates practical experience with:

Bulk RNA-seq
      |
      v
Data validation
      |
      v
Python-based QC
      |
      v
PCA & correlation
      |
      v
DESeq2
      |
      v
Experimental design
      |
      v
GO enrichment
      |
      v
Reactome pathways
      |
      v
Visualization
      |
      v
R Shiny
      |
      v
Git/GitHub reproducibility

The project therefore focuses not only on producing plots, but on connecting:

experimental design → statistical analysis → biological interpretation → interactive exploration

Future Extensions

Potential future extensions include:

Raw FASTQ processing
FastQC / MultiQC integration
Read alignment or pseudoalignment
Transcript quantification
Automated workflow orchestration with Nextflow
Containerized environments
Expanded interactive exploration
Additional public RNA-seq datasets
Independent biological validation

These are outside the current scope of the project.

Author

Rajat Kandpal

B.Tech Biotechnology
Computational Biotechnology

Interested in:

Bioinformatics
Computational Biology
RNA-seq Analysis
Biological Data Science
AI for Healthcare and Life Sciences
Acknowledgements

This project uses publicly available data from the NCBI Gene Expression Omnibus.

The analysis relies on established open-source tools and Bioconductor packages including DESeq2, clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, ggplot2, and Shiny.

License

This repository is intended as an educational and portfolio project.

Please refer to the licenses and citation requirements of the underlying datasets, software packages, and databases before redistributing or using the project for research or commercial purposes.


### Why this version is better

Your actual repository already supports the important technical claims here: the DE script validates counts/metadata and uses `~ donor + treatment`; the QC script performs sample/gene QC; the enrichment and Shiny components are present in the repo. :contentReference[oaicite:2]{index=2}

Also, I intentionally **didn't add fake "key biological discoveries."** Your repository contains the analysis outputs, but the README should not invent scientific conclusions just to look impressive. We can extract the actual strongest findings from the result tables next and add a genuine **Key Findings** section.

And the limitation about not processing FASTQ is important: mature RNA-seq repositories such as `nf-core/rnaseq` explicitly cover FASTQ/BAM input, QC, trimming, alignment/quantification and extensive QC, which is a different scope from what your repository currently implements. :contentReference[oaicite:3]{index=3}

**This README is ready to paste now.** After you paste it, the next thing I'd fix is **not more code**—I'd add 3–4 actual project screenshots/results to the README and then extract your real biological findings from the CSV outputs. That will improve the GitHub presentation more than adding another unnecessary feature.
