from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from modules.module_metabolomics.statistics import bh_fdr


def _is_balanced(metadata: pd.DataFrame, factor_a: str, factor_b: str) -> bool:
    cells = metadata.groupby([factor_a, factor_b]).size()
    expected = metadata[factor_a].nunique() * metadata[factor_b].nunique()
    return len(cells) == expected and cells.nunique() == 1


def _effect_key(index, factor_a, factor_b):
    text = str(index).replace(" ", "")
    if ":" in text:
        return "AxB"
    if factor_a in text:
        return "A"
    if factor_b in text:
        return "B"
    return None


def _fit_feature(frame, factor_a, factor_b, balanced):
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    if balanced:
        formula = f'value ~ C(Q("{factor_a}")) * C(Q("{factor_b}"))'
        model = ols(formula, data=frame).fit()
        table = sm.stats.anova_lm(model, typ=2)
        model_type = "Type-II"
    else:
        formula = f'value ~ C(Q("{factor_a}"), Sum) * C(Q("{factor_b}"), Sum)'
        model = ols(formula, data=frame).fit()
        table = sm.stats.anova_lm(model, typ=3)
        model_type = "Sum contrasts + Type-III"
    return model, table, model_type


def _model_contrast(model, row_a: dict, row_b: dict):
    from patsy import build_design_matrices
    new = pd.DataFrame([row_a, row_b])
    design = build_design_matrices([model.model.data.design_info], new)[0]
    contrast = np.asarray(design[1] - design[0], dtype=float)
    estimate = float(contrast @ model.params)
    covariance = np.asarray(model.cov_params(), dtype=float)
    se = float(np.sqrt(max(contrast @ covariance @ contrast, 0)))
    df = float(model.df_resid)
    t = estimate / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else np.nan
    critical = stats.t.ppf(0.975, df) if df > 0 else np.nan
    return estimate, se, df, t, p, estimate - critical * se, estimate + critical * se


def run_two_factor(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    factor_a: str,
    factor_b: str,
    interaction_fdr: float = 0.05,
    progress=None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = metadata[metadata["Sample type"].astype(str).str.lower().ne("qc")].copy()
    meta = meta[meta[factor_a].fillna("").astype(str).str.strip().ne("") & meta[factor_b].fillna("").astype(str).str.strip().ne("")]
    if meta[factor_a].nunique() < 2 or meta[factor_b].nunique() < 2:
        raise ValueError("Two-factor analysis requires at least two levels for both factors.")
    cells = meta.groupby([factor_a, factor_b]).size()
    if len(cells) != meta[factor_a].nunique() * meta[factor_b].nunique() or (cells < 2).any():
        raise ValueError("Two-factor analysis requires a complete factorial design with at least two replicates per cell.")
    balanced = _is_balanced(meta, factor_a, factor_b)
    sample_ids = [x for x in meta["Sample ID"].astype(str) if x in abundance.columns]
    meta = meta.set_index("Sample ID").loc[sample_ids].reset_index()
    rows = []
    models = {}
    total = len(abundance)
    for number, (feature_id, values) in enumerate(abundance[sample_ids].iterrows(), 1):
        frame = meta.copy(); frame["value"] = pd.to_numeric(values.reindex(sample_ids), errors="coerce").to_numpy()
        frame = frame.dropna(subset=["value"])
        try:
            model, table, model_type = _fit_feature(frame, factor_a, factor_b, balanced)
            residual_ss = float(table.loc["Residual", "sum_sq"])
            effect = {"A": {}, "B": {}, "AxB": {}}
            for idx, row in table.iterrows():
                key = _effect_key(idx, factor_a, factor_b)
                if key:
                    ss = float(row["sum_sq"])
                    effect[key] = {"F": float(row["F"]), "P": float(row["PR(>F)"]), "effect": ss / (ss + residual_ss) if ss + residual_ss > 0 else np.nan}
            rows.append({"Feature ID": str(feature_id), "F_A": effect["A"].get("F"), "P_A": effect["A"].get("P"), "effect_A": effect["A"].get("effect"), "F_B": effect["B"].get("F"), "P_B": effect["B"].get("P"), "effect_B": effect["B"].get("effect"), "F_AxB": effect["AxB"].get("F"), "P_AxB": effect["AxB"].get("P"), "effect_AxB": effect["AxB"].get("effect"), "Model": model_type})
            models[str(feature_id)] = model
        except Exception:
            rows.append({"Feature ID": str(feature_id), "Model": "failed"})
        if callable(progress) and (number == total or number % max(total // 100, 1) == 0): progress(number / total)
    results = pd.DataFrame(rows)
    for effect in ["A", "B", "AxB"]:
        results[f"q_{effect}"] = bh_fdr(results.get(f"P_{effect}", np.nan))

    significant = results.loc[results["q_AxB"] < interaction_fdr, "Feature ID"].astype(str)
    simple_rows = []
    a_levels = list(map(str, pd.unique(meta[factor_a]))); b_levels = list(map(str, pd.unique(meta[factor_b])))
    for feature_id in significant:
        model = models.get(feature_id)
        if model is None: continue
        for b in b_levels:
            for a0, a1 in combinations(a_levels, 2):
                values = _model_contrast(model, {factor_a: a0, factor_b: b}, {factor_a: a1, factor_b: b})
                simple_rows.append({"Feature ID": feature_id, "Effect": "A within B", "Within": b, "Contrast": f"{a1} vs {a0}", "estimate": values[0], "SE": values[1], "df": values[2], "t": values[3], "P": values[4], "CI low": values[5], "CI high": values[6]})
        for a in a_levels:
            for b0, b1 in combinations(b_levels, 2):
                values = _model_contrast(model, {factor_a: a, factor_b: b0}, {factor_a: a, factor_b: b1})
                simple_rows.append({"Feature ID": feature_id, "Effect": "B within A", "Within": a, "Contrast": f"{b1} vs {b0}", "estimate": values[0], "SE": values[1], "df": values[2], "t": values[3], "P": values[4], "CI low": values[5], "CI high": values[6]})
    simple = pd.DataFrame(simple_rows)
    if not simple.empty:
        simple["FDR"] = simple.groupby(["Effect", "Within", "Contrast"])["P"].transform(lambda x: bh_fdr(x.to_numpy()))
    return results, simple

