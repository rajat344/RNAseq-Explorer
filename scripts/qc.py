from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]

COUNTS_FILE = ROOT / "data" / "processed" / "counts.csv"
METADATA_FILE = ROOT / "metadata" / "samples.csv"

RESULTS_DIR = ROOT / "results" / "qc"
FIGURES_DIR = ROOT / "figures" / "qc"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    counts = pd.read_csv(COUNTS_FILE)
    metadata = pd.read_csv(METADATA_FILE)

    sample_columns = [column for column in counts.columns if column != "gene_id"]
    sample_names = [column.removesuffix("_count") for column in sample_columns]

    if counts["gene_id"].duplicated().any():
        raise ValueError("Duplicate gene IDs detected.")

    if not all(
        pd.api.types.is_numeric_dtype(counts[column])
        for column in sample_columns
    ):
        raise ValueError("Count matrix contains non-numeric sample values.")

    if (counts[sample_columns] < 0).any().any():
        raise ValueError("Negative counts detected.")

    if set(sample_names) != set(metadata["sample"]):
        raise ValueError("Count matrix and metadata samples do not match.")

    return counts, metadata, sample_columns, sample_names


def sample_qc(counts, metadata, sample_columns, sample_names):
    sample_counts = counts[sample_columns]

    qc = pd.DataFrame({
        "sample": sample_names,
        "total_counts": sample_counts.sum(axis=0).values,
        "detected_genes": (sample_counts > 0).sum(axis=0).values
    })

    qc = qc.merge(metadata, on="sample", how="left")
    qc.to_csv(RESULTS_DIR / "sample_qc.csv", index=False)

    return qc


def gene_qc(counts, sample_columns):
    sample_counts = counts[sample_columns]

    gene_qc = pd.DataFrame({
        "gene_id": counts["gene_id"],
        "total_count": sample_counts.sum(axis=1),
        "samples_detected": (sample_counts > 0).sum(axis=1),
        "samples_with_10_or_more": (sample_counts >= 10).sum(axis=1)
    })

    gene_qc.to_csv(RESULTS_DIR / "gene_qc.csv", index=False)

    summary = {
        "total_genes": len(gene_qc),
        "zero_count_genes": int((gene_qc["total_count"] == 0).sum()),
        "genes_detected_in_at_least_one_sample": int(
            (gene_qc["samples_detected"] >= 1).sum()
        ),
        "genes_with_10_or_more_counts_in_at_least_two_samples": int(
            (gene_qc["samples_with_10_or_more"] >= 2).sum()
        )
    }

    pd.Series(summary).to_csv(
        RESULTS_DIR / "gene_qc_summary.csv",
        header=["value"]
    )

    return gene_qc, summary


def make_plots(qc):
    ordered = qc.sort_values("total_counts", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.bar(ordered["sample"], ordered["total_counts"])
    plt.xticks(rotation=90)
    plt.ylabel("Total raw counts")
    plt.xlabel("Sample")
    plt.title("RNA-seq library size")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "library_size.png", dpi=300)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.bar(ordered["sample"], ordered["detected_genes"])
    plt.xticks(rotation=90)
    plt.ylabel("Detected genes")
    plt.xlabel("Sample")
    plt.title("Genes detected per sample")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "detected_genes.png", dpi=300)
    plt.close()


def main():
    counts, metadata, sample_columns, sample_names = load_data()

    qc = sample_qc(
        counts,
        metadata,
        sample_columns,
        sample_names
    )

    _, summary = gene_qc(
        counts,
        sample_columns
    )

    make_plots(qc)

    print(f"Genes: {len(counts):,}")
    print(f"Samples: {len(sample_columns)}")
    print(f"Zero-count genes: {summary['zero_count_genes']:,}")
    print(
        "Genes with >=10 counts in >=2 samples: "
        f"{summary['genes_with_10_or_more_counts_in_at_least_two_samples']:,}"
    )
    print(f"QC results: {RESULTS_DIR}")
    print(f"QC figures: {FIGURES_DIR}")


if __name__ == "__main__":
    main()