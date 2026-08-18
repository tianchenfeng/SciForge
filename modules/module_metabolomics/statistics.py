from __future__ import annotations

import numpy as np
import pandas as pd


def bh_fdr(p_values) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    result[valid] = restored
    return result


def holm_adjust(p_values) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values)
    p = values[valid]
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.maximum.accumulate((len(ranked) - np.arange(len(ranked))) * ranked)
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0, 1)
    result[valid] = restored
    return result


def fold_change(mean_reference, mean_target, scale: str) -> tuple[np.ndarray, np.ndarray]:
    ref = np.asarray(mean_reference, dtype=float)
    target = np.asarray(mean_target, dtype=float)
    scale = str(scale).lower()
    if scale == "raw":
        fc = np.divide(target, ref, out=np.full_like(target, np.nan), where=ref != 0)
        log2fc = np.log2(fc, where=fc > 0, out=np.full_like(fc, np.nan))
        return fc, log2fc
    difference = target - ref
    if scale == "log2":
        return np.power(2.0, difference), difference
    if scale == "log10":
        fc = np.power(10.0, difference)
        return fc, np.log2(fc)
    return np.full_like(target, np.nan), np.full_like(target, np.nan)


def fold_change_rule(scale: str, threshold: float, pseudocount: float | None = None) -> dict:
    scale = str(scale).lower()
    common = {
        "computed": scale in {"raw", "log2", "log10"},
        "input_scale": scale,
        "contrast_direction": "target_over_reference",
        "contrast_label_rule": "TARGET vs REFERENCE",
        "interpretation": "FC > 1 means TARGET abundance is higher than REFERENCE abundance.",
        "pseudocount": pseudocount,
        "threshold": {"upper": float(threshold), "lower": float(1.0 / threshold)},
    }
    if scale == "raw":
        common.update({
            "fc_formula": "mean(target) / mean(reference)",
            "log2fc_formula": "log2(mean(target) / mean(reference))",
            "zero_handling": "FC and log2FC are unavailable when mean(reference) is zero; no pseudocount is added for FC.",
        })
    elif scale == "log2":
        common.update({
            "fc_formula": "2 ** (mean_log2(target) - mean_log2(reference))",
            "log2fc_formula": "mean_log2(target) - mean_log2(reference)",
            "zero_handling": "Not applicable on the confirmed log2 scale.",
        })
    elif scale == "log10":
        common.update({
            "fc_formula": "10 ** (mean_log10(target) - mean_log10(reference))",
            "log2fc_formula": "log2(10 ** (mean_log10(target) - mean_log10(reference)))",
            "zero_handling": "Not applicable on the confirmed log10 scale.",
        })
    else:
        common.update({
            "reason": "Fold change is disabled because the abundance scale is unknown.",
            "fc_formula": None,
            "log2fc_formula": None,
            "zero_handling": None,
        })
    return common


def summarize_significant_features(differential: pd.DataFrame) -> dict:
    if differential.empty or "Significance class" not in differential:
        return {"records": 0, "unique_features": 0, "up": 0, "down": 0, "per_contrast": {}}
    significant = differential[differential["Significance class"].isin(["Up", "Down"])].copy()
    if significant.empty:
        return {"records": 0, "unique_features": 0, "up": 0, "down": 0, "per_contrast": {}}
    feature_id = significant["Feature ID"].fillna("").astype(str)
    if "Feature Key" in significant:
        stable_key = significant["Feature Key"].fillna("").astype(str)
        stable_key = stable_key.where(stable_key.ne(""), feature_id)
    elif "Ion mode" in significant:
        stable_key = significant["Ion mode"].fillna("unknown").astype(str) + "::" + feature_id
    else:
        stable_key = feature_id
    per_contrast = significant.groupby("Contrast").size().astype(int).to_dict() if "Contrast" in significant else {}
    return {
        "records": int(len(significant)),
        "unique_features": int(stable_key.nunique()),
        "up": int(significant["Significance class"].eq("Up").sum()),
        "down": int(significant["Significance class"].eq("Down").sum()),
        "per_contrast": per_contrast,
    }


def significance_class(fdr, fc, fdr_threshold=0.05, fc_threshold=1.2) -> np.ndarray:
    fdr = np.asarray(fdr, dtype=float)
    fc = np.asarray(fc, dtype=float)
    result = np.full(fdr.shape, "NS", dtype=object)
    result[(fdr < fdr_threshold) & (fc >= fc_threshold)] = "Up"
    result[(fdr < fdr_threshold) & (fc <= 1.0 / fc_threshold)] = "Down"
    result[(fdr < fdr_threshold) & ~np.isfinite(fc)] = "Statistically significant; FC unavailable"
    return result
