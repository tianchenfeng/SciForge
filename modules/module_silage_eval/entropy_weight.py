from __future__ import annotations

import numpy as np
import pandas as pd


def entropy_weights(normalized_df: pd.DataFrame) -> pd.Series:
    data = normalized_df.astype(float).clip(lower=0.0)
    if data.empty:
        return pd.Series(dtype=float)

    variable_mask = data.max(axis=0) - data.min(axis=0) > 1e-12
    col_sums = data.sum(axis=0).replace(0, np.nan)
    p = data.div(col_sums, axis=1).fillna(0.0)
    n = len(data)
    if n <= 1:
        return pd.Series(1.0 / len(data.columns), index=data.columns)
    k = 1.0 / np.log(n)
    entropy = -k * (p * np.log(p.replace(0, np.nan))).sum(axis=0).fillna(0.0)
    diversity = 1.0 - entropy
    diversity = diversity.where(variable_mask, 0.0)
    if float(diversity.sum()) <= 0:
        return pd.Series(1.0 / len(data.columns), index=data.columns)
    return diversity / diversity.sum()
