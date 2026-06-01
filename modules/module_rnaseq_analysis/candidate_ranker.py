from __future__ import annotations

import numpy as np
import pandas as pd

from modules.module_rnaseq_analysis.gene_set_ops import keyword_hit_count


def rank_candidates(
    deg_with_annotation: pd.DataFrame,
    pathway_gene_map: pd.DataFrame | None = None,
    keywords: list[str] | None = None,
) -> pd.DataFrame:
    result = deg_with_annotation.copy()
    if result.empty:
        result["KeywordHits"] = []
        result["PathwayHits"] = []
        result["CandidateScore"] = []
        return result

    keywords = keywords or []
    pathway_gene_map = pathway_gene_map if pathway_gene_map is not None else pd.DataFrame()
    pathway_counts = (
        pathway_gene_map.groupby("Gene").size().rename("PathwayHits")
        if not pathway_gene_map.empty and "Gene" in pathway_gene_map.columns
        else pd.Series(dtype=float, name="PathwayHits")
    )

    result["KeywordHits"] = result.apply(lambda row: keyword_hit_count(row, keywords), axis=1)
    result = result.merge(pathway_counts, left_on="Gene", right_index=True, how="left")
    result["PathwayHits"] = result["PathwayHits"].fillna(0).astype(int)
    padj = result["padj"].clip(lower=1e-300)
    result["CandidateScore"] = (
        result["log2FC"].abs() * 10.0
        + (-np.log10(padj))
        + result["KeywordHits"] * 5.0
        + result["PathwayHits"] * 1.5
    )
    return result.sort_values(
        ["KeywordHits", "CandidateScore", "padj", "log2FC"],
        ascending=[False, False, True, False],
    ).reset_index(drop=True)
