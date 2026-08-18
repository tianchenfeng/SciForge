from __future__ import annotations

import numpy as np
import pandas as pd

from modules.module_metabolomics.data_model import MetabolomicsDataset


def compute_qc(dataset: MetabolomicsDataset) -> dict[str, pd.DataFrame | dict]:
    values = dataset.abundance.apply(pd.to_numeric, errors="coerce")
    missing_by_sample = values.isna().mean(axis=0).rename("Missing fraction").rename_axis("Sample ID").reset_index()
    missing_by_feature = values.isna().mean(axis=1).rename("Missing fraction").rename_axis("Feature ID").reset_index()
    sample_summary = pd.DataFrame({
        "Sample ID": values.columns,
        "Median": values.median(axis=0).to_numpy(),
        "Mean": values.mean(axis=0).to_numpy(),
        "Missing fraction": values.isna().mean(axis=0).to_numpy(),
    }).merge(dataset.samples, on="Sample ID", how="left")

    qc_ids = [sample for sample in dataset.qc_samples if sample in values.columns]
    if not qc_ids:
        return {
            "missing_by_sample": missing_by_sample,
            "missing_by_feature": missing_by_feature,
            "sample_summary": sample_summary,
            "qc_table": pd.DataFrame(),
            "qc_correlation": pd.DataFrame(),
            "summary": {"QC samples": 0, "Message": "No QC samples were provided; QC-based stability assessment is unavailable."},
        }

    qc = values[qc_ids]
    mean = qc.mean(axis=1); sd = qc.std(axis=1, ddof=1)
    rsd = sd / mean.abs().replace(0, np.nan)
    feature_info = dataset.features.set_index("Feature ID") if "Feature ID" in dataset.features else pd.DataFrame()
    ion = feature_info.get("Ion mode", pd.Series(index=values.index, dtype=object)).reindex(values.index)
    qc_table = pd.DataFrame({"Feature ID": values.index.astype(str), "Ion mode": ion.to_numpy(), "QC mean": mean.to_numpy(), "QC SD": sd.to_numpy(), "QC RSD": rsd.to_numpy(), "Pass20": (rsd < 0.2).to_numpy(), "Pass30": (rsd < 0.3).to_numpy()})
    correlation = qc.corr(method="pearson")
    upper = correlation.where(np.triu(np.ones(correlation.shape), 1).astype(bool)).stack()
    summary = {
        "QC samples": len(qc_ids),
        "RSD <20%": float((rsd < 0.2).mean()),
        "RSD <30%": float((rsd < 0.3).mean()),
        "QC mean correlation": float(upper.mean()) if len(upper) else np.nan,
    }
    return {"missing_by_sample": missing_by_sample, "missing_by_feature": missing_by_feature, "sample_summary": sample_summary, "qc_table": qc_table, "qc_correlation": correlation, "summary": summary}
