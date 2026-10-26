from __future__ import annotations

import numpy as np
import pandas as pd

from modules.module_metabolomics.data_model import MetabolomicsDataset


def _row(check: str, status: str, value, message: str) -> dict:
    return {"Check": check, "Status": status, "Value": value, "Message": message}


def validate_dataset(dataset: MetabolomicsDataset) -> pd.DataFrame:
    abundance = dataset.abundance
    samples = dataset.samples.copy()
    features = dataset.features.copy()
    rows: list[dict] = []

    rows.append(_row("Feature count", "PASS" if len(abundance) else "ERROR", len(abundance), "Features in the analysis matrix."))
    rows.append(_row("Sample count", "PASS" if abundance.shape[1] else "ERROR", abundance.shape[1], "Samples in the analysis matrix."))
    rows.append(_row("Biological samples", "PASS" if dataset.biological_samples else "ERROR", len(dataset.biological_samples), "Non-QC samples."))
    rows.append(_row("QC samples", "PASS" if dataset.qc_samples else "WARNING", len(dataset.qc_samples), "QC-based stability metrics require QC samples."))

    sample_duplicates = int(samples["Sample ID"].astype(str).duplicated().sum())
    feature_duplicates = int(abundance.index.astype(str).duplicated().sum())
    rows.append(_row("Duplicated Sample ID", "ERROR" if sample_duplicates else "PASS", sample_duplicates, "Sample IDs must be unique."))
    rows.append(_row("Duplicated feature ID", "ERROR" if feature_duplicates else "PASS", feature_duplicates, "Feature IDs must be unique within the imported matrix."))

    known_samples = set(samples["Sample ID"].astype(str))
    missing_metadata = sorted(set(map(str, abundance.columns)) - known_samples)
    rows.append(_row("Sample metadata alignment", "ERROR" if missing_metadata else "PASS", len(missing_metadata), "All abundance columns must occur in sample metadata."))

    bio = samples[samples["Sample type"].astype(str).str.lower().ne("qc")]
    empty_group = int(bio["Group"].fillna("").astype(str).str.strip().eq("").sum())
    rows.append(_row("Empty biological group", "ERROR" if empty_group else "PASS", empty_group, "Confirm grouping before formal analysis."))
    group_counts = bio["Group"].astype(str).value_counts().to_dict()
    rows.append(_row("Replicates per group", "WARNING" if any(v < 3 for v in group_counts.values()) else "PASS", str(group_counts), "At least three replicates per group are recommended."))

    numeric = abundance.apply(pd.to_numeric, errors="coerce")
    magnitudes = np.abs(numeric.to_numpy(dtype=float))
    infinite_count = int(np.isinf(magnitudes).sum())
    unsafe_count = int(((magnitudes > 1e150) | ((magnitudes > 0) & (magnitudes < 1e-150))).sum())
    rows.append(_row("Finite abundance", "ERROR" if infinite_count else "PASS", infinite_count,
                     "Infinity is not a valid abundance; check the source data."))
    rows.append(_row("Numeric range", "ERROR" if unsafe_count else "PASS", unsafe_count,
                     "Nonzero abundance must be between 1e-150 and 1e150; convert units before analysis."))
    total = max(numeric.size, 1)
    missing_rate = float(numeric.isna().to_numpy().sum() / total)
    negative_rate = float((numeric < 0).to_numpy().sum() / total)
    zero_rate = float((numeric == 0).to_numpy().sum() / total)
    rows.append(_row("Missing values", "WARNING" if missing_rate else "PASS", missing_rate, "Fraction of missing abundance cells."))
    rows.append(_row("Negative values", "WARNING" if negative_rate else "PASS", negative_rate, "Negative values affect logarithmic transformations."))
    rows.append(_row("Zero proportion", "WARNING" if zero_rate > 0.2 else "PASS", zero_rate, "Zeros require an explicit missing/transformation policy."))

    if "Ion mode" in features:
        modes = sorted(set(features["Ion mode"].dropna().astype(str)) - {"", "-", "unknown"})
    else:
        modes = []
    rows.append(_row("Ion mode", "PASS" if modes else "WARNING", ", ".join(modes) or "unknown", "Detected ionization modes."))

    key_dups = 0
    if {"Feature ID", "Ion mode"}.issubset(features.columns):
        key_dups = int(features[["Ion mode", "Feature ID"]].astype(str).duplicated().sum())
    rows.append(_row("Duplicated ion_mode + feature_id", "ERROR" if key_dups else "PASS", key_dups, "Composite feature key must be unique."))

    if "Identified" in features:
        annotated = int(features["Identified"].fillna(False).astype(bool).sum())
    elif "Annotation level" in features:
        annotated = int((~features["Annotation level"].fillna("").astype(str).isin(["", "-", "unknown"])).sum())
    else:
        annotated = 0
    kegg = int((~features.get("KEGG", pd.Series(index=features.index, dtype=object)).fillna("").astype(str).isin(["", "-", "nan"])).sum())
    rows.append(_row("Annotation mapping rate", "WARNING" if not annotated else "PASS", annotated / max(len(features), 1), f"{annotated} annotated features."))
    rows.append(_row("KEGG mapping rate", "WARNING" if not kegg else "PASS", kegg / max(len(features), 1), f"{kegg} features with explicit KEGG IDs."))
    rows.append(_row("Abundance scale", "WARNING" if dataset.scale == "unknown" else "PASS", dataset.scale, "Unknown scale must be confirmed before fold-change calculation."))
    return pd.DataFrame(rows)


def has_blocking_errors(validation: pd.DataFrame) -> bool:
    return validation["Status"].astype(str).eq("ERROR").any()
