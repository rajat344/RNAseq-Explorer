from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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

    return np.log2(cpm + 1)


def calculate_correlation(expression):
    return expression.corr(method="pearson")


def save_correlation_matrix(correlation):
    correlation.to_csv(
        RESULTS_DIR / "sample_correlation.csv"
    )


def plot_heatmap(correlation, metadata):
    sample_order = metadata["sample"].tolist()

    correlation = correlation.loc[
        sample_order,
        sample_order
    ]

    plt.figure(figsize=(11, 9))

    plt.imshow(
        correlation,
        aspect="auto",
        interpolation="nearest"
    )

    plt.colorbar(label="Pearson correlation")

    plt.xticks(
        range(len(sample_order)),
        sample_order,
        rotation=90
    )

    plt.yticks(
        range(len(sample_order)),
        sample_order
    )

    plt.xlabel("Sample")
    plt.ylabel("Sample")
    plt.title("Sample-to-sample RNA-seq correlation")

    plt.tight_layout()

    plt.savefig(
        FIGURES_DIR / "sample_correlation_heatmap.png",
        dpi=300
    )

    plt.close()


def main():
    expression, metadata = load_data()

    log_cpm = normalize_counts(expression)

    correlation = calculate_correlation(log_cpm)

    save_correlation_matrix(correlation)

    plot_heatmap(correlation, metadata)

    off_diagonal = correlation.values[
        ~np.eye(len(correlation), dtype=bool)
    ]

    print(f"Mean pairwise correlation: {off_diagonal.mean():.3f}")
    print(f"Minimum pairwise correlation: {off_diagonal.min():.3f}")
    print(f"Maximum pairwise correlation: {off_diagonal.max():.3f}")

    print(
        f"Correlation matrix: "
        f"{RESULTS_DIR / 'sample_correlation.csv'}"
    )

    print(
        f"Correlation heatmap: "
        f"{FIGURES_DIR / 'sample_correlation_heatmap.png'}"
    )


if __name__ == "__main__":
    main()