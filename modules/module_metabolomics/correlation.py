from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from modules.module_metabolomics.adapters.base import read_table_smart
from modules.module_metabolomics.statistics import bh_fdr


def load_phenotype(path) -> pd.DataFrame:
    table = read_table_smart(path, low_memory=False)
    lower = {str(c).strip().lower(): c for c in table.columns}
    sample_col = lower.get("sample id") or lower.get("sample_id") or lower.get("sample")
    if sample_col is None:
        raise ValueError("Phenotype table requires a Sample ID column.")
    return table.rename(columns={sample_col: "Sample ID"})


def correlate_features(
    abundance: pd.DataFrame,
    phenotype: pd.DataFrame,
    method: str = "spearman",
    feature_ids: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    phenotype = phenotype.copy(); phenotype["Sample ID"] = phenotype["Sample ID"].astype(str)
    samples = [s for s in phenotype["Sample ID"] if s in abundance.columns]
    if len(samples) < 3:
        raise ValueError("Correlation requires at least three matched Sample IDs.")
    phenotype = phenotype.set_index("Sample ID").loc[samples]
    phenotypes = [c for c in phenotype if pd.to_numeric(phenotype[c], errors="coerce").notna().sum() >= 3]
    if not phenotypes:
        raise ValueError("No numeric phenotype columns with at least three observations were found.")
    matrix = abundance.loc[feature_ids] if feature_ids else abundance
    matrix = matrix[samples].apply(pd.to_numeric, errors="coerce")
    rows = []
    for feature_id, values in matrix.iterrows():
        x = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
        for phenotype_name in phenotypes:
            y = pd.to_numeric(phenotype[phenotype_name], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 3 or np.nanstd(x[valid]) == 0 or np.nanstd(y[valid]) == 0:
                coefficient, p = np.nan, np.nan
            elif method.lower() == "pearson":
                coefficient, p = stats.pearsonr(x[valid], y[valid])
            else:
                coefficient, p = stats.spearmanr(x[valid], y[valid])
            rows.append({"Feature ID": str(feature_id), "Phenotype": phenotype_name, "Correlation": coefficient, "P": p, "N": int(valid.sum()), "Method": method})
    long = pd.DataFrame(rows); long["FDR"] = bh_fdr(long["P"])
    corr = long.pivot(index="Feature ID", columns="Phenotype", values="Correlation")
    q = long.pivot(index="Feature ID", columns="Phenotype", values="FDR")
    return long.sort_values(["FDR", "P"]), corr, q


def correlate_16s_metabolomics(taxon_abundance: pd.DataFrame, metabolite_abundance: pd.DataFrame) -> pd.DataFrame:
    common = [sample for sample in taxon_abundance.columns if sample in metabolite_abundance.columns]
    if len(common) < 3:
        raise ValueError("16S × metabolomics correlation requires at least three matched Sample IDs.")
    rows = []
    for taxon, tx in taxon_abundance[common].iterrows():
        x = pd.to_numeric(tx, errors="coerce").to_numpy(dtype=float)
        for feature, mt in metabolite_abundance[common].iterrows():
            y = pd.to_numeric(mt, errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 3 or np.std(x[valid]) == 0 or np.std(y[valid]) == 0: coefficient, p = np.nan, np.nan
            else: coefficient, p = stats.spearmanr(x[valid], y[valid])
            rows.append({"Taxon": str(taxon), "Feature ID": str(feature), "Spearman": coefficient, "P": p, "N": int(valid.sum())})
    result = pd.DataFrame(rows); result["FDR"] = bh_fdr(result["P"])
    return result.sort_values(["FDR", "P"])

