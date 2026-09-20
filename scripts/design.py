from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]

METADATA_FILE = ROOT / "metadata" / "samples.csv"
RESULTS_DIR = ROOT / "results" / "design"
FIGURES_DIR = ROOT / "figures" / "design"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_metadata():
    metadata = pd.read_csv(METADATA_FILE)

    required_columns = {
        "sample",
        "geo_accession",
        "donor",
        "microenvironment",
        "treatment",
    }

    missing_columns = required_columns.difference(metadata.columns)

    if missing_columns:
        names = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing metadata columns: {names}")

    if metadata["sample"].duplicated().any():
        raise ValueError("Duplicate sample names detected.")

    if metadata[["donor", "microenvironment", "treatment"]].isna().any().any():
        raise ValueError("Missing values found in experimental metadata.")

    return metadata


def create_group_summary(metadata):
    summary = (
        metadata
        .groupby(
            ["microenvironment", "treatment"],
            sort=True
        )
        .size()
        .reset_index(name="sample_count")
    )

    summary.to_csv(
        RESULTS_DIR / "group_summary.csv",
        index=False
    )

    return summary


def create_donor_summary(metadata):
    donor_summary = (
        metadata
        .groupby(
            ["donor", "microenvironment", "treatment"],
            sort=True
        )
        .size()
        .reset_index(name="sample_count")
    )

    donor_summary.to_csv(
        RESULTS_DIR / "donor_summary.csv",
        index=False
    )

    return donor_summary


def create_design_matrix(metadata):
    design = pd.get_dummies(
        metadata[
            [
                "microenvironment",
                "treatment",
                "donor"
            ]
        ],
        drop_first=False,
        dtype=int
    )

    design.insert(
        0,
        "sample",
        metadata["sample"].values
    )

    design.to_csv(
        RESULTS_DIR / "design_matrix.csv",
        index=False
    )

    return design


def create_contrast_plan():
    comparisons = [
        {
            "microenvironment": "Microenvironment_1",
            "treatment": "5-FU",
            "reference": "Control",
        },
        {
            "microenvironment": "Microenvironment_1",
            "treatment": "Oxaliplatin",
            "reference": "Control",
        },
        {
            "microenvironment": "Microenvironment_2",
            "treatment": "5-FU",
            "reference": "Control",
        },
        {
            "microenvironment": "Microenvironment_2",
            "treatment": "Oxaliplatin",
            "reference": "Control",
        },
        {
            "microenvironment": "Microenvironment_3",
            "treatment": "Immunotherapy",
            "reference": "Control",
        },
        {
            "microenvironment": "Microenvironment_3",
            "treatment": "Chemotherapy + Immunotherapy",
            "reference": "Control",
        },
    ]

    contrasts = pd.DataFrame(comparisons)

    contrasts["comparison"] = (
        contrasts["treatment"]
        + " vs "
        + contrasts["reference"]
    )

    contrasts = contrasts[
        [
            "microenvironment",
            "comparison",
            "treatment",
            "reference"
        ]
    ]

    contrasts.to_csv(
        RESULTS_DIR / "contrast_plan.csv",
        index=False
    )

    return contrasts


def create_design_figure(summary):
    table = summary.pivot(
        index="microenvironment",
        columns="treatment",
        values="sample_count"
    ).fillna(0)

    ax = table.plot(
        kind="bar",
        figsize=(10, 6)
    )

    ax.set_xlabel("Microenvironment")
    ax.set_ylabel("Number of samples")
    ax.set_title("RNA-seq experimental design")
    ax.legend(title="Treatment")

    plt.tight_layout()

    plt.savefig(
        FIGURES_DIR / "sample_design.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


def main():
    metadata = load_metadata()

    summary = create_group_summary(metadata)
    create_donor_summary(metadata)
    create_design_matrix(metadata)
    contrasts = create_contrast_plan()

    create_design_figure(summary)

    print(f"Samples: {len(metadata)}")
    print(f"Donors: {metadata['donor'].nunique()}")
    print(f"Microenvironments: {metadata['microenvironment'].nunique()}")
    print(f"Treatments: {metadata['treatment'].nunique()}")

    print("\nExperimental design:")
    print(summary.to_string(index=False))

    print("\nPlanned contrasts:")
    print(contrasts.to_string(index=False))

    print(f"\nDesign results: {RESULTS_DIR}")
    print(f"Design figure: {FIGURES_DIR / 'sample_design.png'}")


if __name__ == "__main__":
    main()