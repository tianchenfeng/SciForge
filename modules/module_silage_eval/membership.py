from __future__ import annotations

import pandas as pd


NEGATIVE_KEYWORDS = ["pH", "NH3-N", "AN/TN", "BA", "DML", "NDF", "ADF"]
POSITIVE_KEYWORDS = ["LA", "AA", "PA", "CP", "WSC", "DM", "IVDMD", "NDFD", "ADFD", "MCP"]


def _direction_patterns() -> list[tuple[str, str]]:
    patterns = [(key.upper(), "negative") for key in NEGATIVE_KEYWORDS]
    patterns.extend((key.upper(), "positive") for key in POSITIVE_KEYWORDS)
    return sorted(patterns, key=lambda item: len(item[0]), reverse=True)


def infer_indicator_directions(columns) -> dict[str, str]:
    directions = {}
    for col in columns:
        text = str(col).upper()
        direction = "positive"
        for key, matched_direction in _direction_patterns():
            if key in text:
                direction = matched_direction
                break
        directions[str(col)] = direction
    return directions


def range_normalize(df: pd.DataFrame, directions: dict[str, str]) -> pd.DataFrame:
    normalized = pd.DataFrame(index=df.index)
    for col in df.columns:
        values = pd.to_numeric(df[col], errors="coerce").astype(float)
        min_v = values.min()
        max_v = values.max()
        span = max_v - min_v
        if pd.isna(span) or span == 0:
            normalized[col] = 0.0
        elif directions.get(str(col), "positive") == "negative":
            normalized[col] = (max_v - values) / span
        else:
            normalized[col] = (values - min_v) / span
    return normalized.fillna(0.0)


def membership_scores(normalized_df: pd.DataFrame, weights: pd.Series) -> pd.Series:
    aligned_weights = weights.reindex(normalized_df.columns).fillna(0.0)
    return normalized_df.mul(aligned_weights, axis=1).sum(axis=1)
