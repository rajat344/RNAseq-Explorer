from pathlib import Path
import gzip
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

MATRIX_FILE = ROOT / "data" / "raw" / "GSE335921_raw_gene_matrix.txt.gz"
SOFT_FILE = ROOT / "data" / "raw" / "GSE335921_family.soft.gz"

COUNTS_FILE = ROOT / "data" / "processed" / "counts.csv"
METADATA_FILE = ROOT / "metadata" / "samples.csv"


def load_counts():
    data = pd.read_csv(
        MATRIX_FILE,
        sep="\t",
        compression="gzip",
        usecols=lambda column: column == "gene_id" or column.endswith("_count")
    )

    count_columns = [column for column in data.columns if column != "gene_id"]

    if len(count_columns) != 36:
        raise ValueError(f"Expected 36 count columns, found {len(count_columns)}")

    data.to_csv(COUNTS_FILE, index=False)
    return data, count_columns


def load_metadata():
    samples = []
    current = None

    with gzip.open(SOFT_FILE, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip()

            if line.startswith("!Sample_title = "):
                if current and "geo_accession" in current:
                    samples.append(current)
                current = {"title": line.split(" = ", 1)[1]}

            elif current and line.startswith("!Sample_geo_accession = "):
                current["geo_accession"] = line.split(" = ", 1)[1]

            elif current and line.startswith("!Sample_characteristics_ch1 = treatment: "):
                current["treatment"] = line.split(": ", 1)[1]

    if current and "geo_accession" in current:
        samples.append(current)

    if len(samples) != 36:
        raise ValueError(f"Expected 36 samples, found {len(samples)}")

    records = []

    replacements = {
        "CPcontrol": "CPc",
        "CPPD-1": "CPp1",
        "CPPD1": "CPp1",
        "CPSOXPD1": "CPs1",
        "control": "con",
        "5fu": "fu"
    }

    treatment_names = {
        "vehicle": "Control",
        "Fluorouracil": "5-FU",
        "Oxaliplatin": "Oxaliplatin",
        "Immunotherapy": "Immunotherapy",
        "Chemotherapy combined with immunotherapy": "Chemotherapy + Immunotherapy"
    }

    for sample in samples:
        sample_id = sample["title"].removesuffix(" sample").replace(" ", "")

        for old, new in replacements.items():
            sample_id = sample_id.replace(old, new)

        donor = re.match(r"^(P\d+)", sample_id).group(1)

        if "_CP" in sample_id:
            microenvironment = "Microenvironment_3"
        elif "_C" in sample_id:
            microenvironment = "Microenvironment_2"
        else:
            microenvironment = "Microenvironment_1"

        records.append({
            "sample": sample_id,
            "geo_accession": sample["geo_accession"],
            "donor": donor,
            "microenvironment": microenvironment,
            "treatment": treatment_names[sample["treatment"]]
        })

    return pd.DataFrame(records)


counts, count_columns = load_counts()
metadata = load_metadata()

count_samples = {column.removesuffix("_count") for column in count_columns}
metadata_samples = set(metadata["sample"])

if count_samples != metadata_samples:
    missing = count_samples - metadata_samples
    extra = metadata_samples - count_samples
    raise ValueError(f"Sample mismatch. Missing: {missing}. Extra: {extra}")

metadata.to_csv(METADATA_FILE, index=False)

print(f"Genes: {len(counts):,}")
print(f"Samples: {len(metadata)}")
print(f"Count columns validated: {len(count_samples)}")
print(f"Metadata saved: {METADATA_FILE}")
print(f"Counts saved: {COUNTS_FILE}")