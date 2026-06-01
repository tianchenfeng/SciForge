from __future__ import annotations

import pandas as pd


def filter_deg(df: pd.DataFrame, log2fc_threshold: float = 1.0, padj_threshold: float = 0.05) -> pd.DataFrame:
    mask = (df["log2FC"].abs() >= float(log2fc_threshold)) & (df["padj"] <= float(padj_threshold))
    result = df.loc[mask].copy()
    result["Regulation"] = result["log2FC"].apply(lambda value: "Up" if value > 0 else "Down")
    return result.reset_index(drop=True)


def split_up_down(deg_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    if deg_df.empty:
        return deg_df.copy(), deg_df.copy()
    up = deg_df[deg_df["Regulation"] == "Up"].sort_values(["padj", "log2FC"], ascending=[True, False])
    down = deg_df[deg_df["Regulation"] == "Down"].sort_values(["padj", "log2FC"], ascending=[True, True])
    return up.reset_index(drop=True), down.reset_index(drop=True)
