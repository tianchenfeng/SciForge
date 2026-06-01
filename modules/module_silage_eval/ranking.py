from __future__ import annotations

import pandas as pd


def build_ranking(treatments: pd.Series, scores: pd.Series) -> pd.DataFrame:
    result = pd.DataFrame({
        "Treatment": treatments.astype(str).values,
        "Score": scores.astype(float).values,
        "_order": range(len(treatments)),
    })
    result = result.sort_values(["Score", "_order"], ascending=[False, True], kind="mergesort").reset_index(drop=True)
    result["Rank"] = range(1, len(result) + 1)
    return result.drop(columns=["_order"])
