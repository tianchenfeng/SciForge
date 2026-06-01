from __future__ import annotations

import pandas as pd


def merge_annotation(deg_df: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    if annotation.empty:
        return deg_df.copy()
    return deg_df.merge(annotation, on="Gene", how="left")
