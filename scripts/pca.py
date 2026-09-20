from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA


ROOT = Path(__file__).resolve().parents[1]

COUNTS_FILE = ROOT / "data" / "processed" / "counts.csv"
METADATA_FILE = ROOT / "metadata" / "samples.csv"

RESULTS_DIR = ROOT / "results" / "qc"
FIGURES_DIR = ROOT / "figures" / "qc"


def load_data():
    counts = pd.read_csv(COUNTS_FILE)
    metadata = pd.read_csv(METADATA_FILE)

    sample_columns = [
        column for column in counts.columns
        if column != "gene_id"
    ]

    sample_names = [
        column.removesuffix("_count")
        for column in sample_columns
    ]

    expression = counts[sample_columns].copy()
    expression.columns = sample_names

    return expression, metadata


def normalize_counts(expression):
    library_sizes = expression.sum(axis=0)

    cpm = expression.div(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 1)

    return log_cpm


def run_pca(log_cpm):
    gene_variance = log_cpm.var(axis=1)

    variable_genes = gene_variance.nlargest(
        min(5000, len(gene_variance))
    ).index

    selected = log_cpm.loc[variable_genes].T

    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(selected)

    result = pd.DataFrame(
        coordinates,
        index=selected.index,
        columns=["PC1", "PC2"]
    )

    explained = pca.explained_variance_ratio_

    return result, explained


def plot_pca(pca_result, explained):
    plt.figure(figsize=(10, 7))

    for microenvironment in pca_result["microenvironment"].unique():
        subset = pca_result[
            pca_result["microenvironment"] == microenvironment
        ]

        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            label=microenvironment,
            s=65
        )

    pc1 = explained[0] * 100
    pc2 = explained[1] * 100

    plt.xlabel(f"PC1 ({pc1:.1f}% variance)")
    plt.ylabel(f"PC2 ({pc2:.1f}% variance)")
    plt.title("PCA of RNA-seq samples")
    plt.legend(title="Microenvironment")
    plt.tight_layout()

    plt.savefig(
        FIGURES_DIR / "pca_samples.png",
        dpi=300
    )

    plt.close()


def main():
    expression, metadata = load_data()

    log_cpm = normalize_counts(expression)

    pca_result, explained = run_pca(log_cpm)

    pca_result = pca_result.reset_index()
    pca_result = pca_result.rename(columns={"index": "sample"})

    pca_result = pca_result.merge(
        metadata,
        on="sample",
        how="left"
    )

    pca_result.to_csv(
        RESULTS_DIR / "pca_coordinates.csv",
        index=False
    )

    plot_pca(pca_result, explained)

    pc1 = explained[0] * 100
    pc2 = explained[1] * 100

    print(f"PC1 variance explained: {pc1:.2f}%")
    print(f"PC2 variance explained: {pc2:.2f}%")
    print(f"Total variance explained: {pc1 + pc2:.2f}%")
    print(f"PCA results: {RESULTS_DIR / 'pca_coordinates.csv'}")
    print(f"PCA figure: {FIGURES_DIR / 'pca_samples.png'}")


if __name__ == "__main__":
    main()