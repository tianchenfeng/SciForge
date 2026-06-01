from __future__ import annotations

import pandas as pd

from modules.module_rnaseq_analysis.gene_set_ops import keyword_hit_count, split_gene_list
from modules.module_rnaseq_analysis.preprocess import read_table


def standardize_enrichment_table(input_file: str | None, source: str, keywords: list[str], padj_threshold: float) -> dict[str, pd.DataFrame]:
    if not input_file:
        empty = _empty_enrichment(source)
        return {"all": empty, "top": empty, "significant": empty, "keyword": empty}

    df = read_table(input_file).copy()
    lower_map = {str(col).lower(): col for col in df.columns}
    id_col = _first_existing(lower_map, ["id", "go", "kegg", "termid", "category"])
    term_col = _first_existing(lower_map, ["description", "term", "name", "pathway"])
    pvalue_col = _first_existing(lower_map, ["pvalue", "p.value", "p_value", "p"])
    padj_col = _first_existing(lower_map, ["padj", "fdr", "qvalue", "q_value", "p.adjust"])
    gene_col = _first_existing(lower_map, ["geneid", "genes", "gene", "core_enrichment"])

    result = pd.DataFrame({
        "Source": source,
        "ID": df[id_col].astype(str) if id_col else "",
        "Term": df[term_col].astype(str) if term_col else df.iloc[:, 0].astype(str),
        "pvalue": pd.to_numeric(df[pvalue_col], errors="coerce") if pvalue_col else pd.NA,
        "padj": pd.to_numeric(df[padj_col], errors="coerce") if padj_col else pd.NA,
        "Genes": df[gene_col].astype(str) if gene_col else "",
    })
    result["GeneCount"] = result["Genes"].apply(lambda value: len(split_gene_list(value)))
    result["KeywordHits"] = result.apply(lambda row: keyword_hit_count(row, keywords), axis=1)
    result = result.sort_values(["padj", "pvalue"], na_position="last").reset_index(drop=True)
    return {
        "all": result,
        "top": result.head(20).copy(),
        "significant": result[result["padj"] <= float(padj_threshold)].copy(),
        "keyword": result[result["KeywordHits"] > 0].copy(),
    }


def _first_existing(lower_map: dict[str, str], names: list[str]) -> str | None:
    for name in names:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None


def _empty_enrichment(source: str) -> pd.DataFrame:
    return pd.DataFrame(columns=["Source", "ID", "Term", "pvalue", "padj", "Genes", "GeneCount", "KeywordHits"]).assign(Source=source)
