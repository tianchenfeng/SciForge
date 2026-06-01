from __future__ import annotations

import pandas as pd

from modules.module_rnaseq_analysis.gene_set_ops import split_gene_list


def build_pathway_gene_map(go_df: pd.DataFrame, kegg_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, df in [("GO", go_df), ("KEGG", kegg_df)]:
        for _, row in df.iterrows():
            for gene in split_gene_list(row.get("Genes", "")):
                rows.append({
                    "Source": source,
                    "PathwayID": row.get("ID", ""),
                    "Pathway": row.get("Term", ""),
                    "Gene": gene,
                    "padj": row.get("padj", pd.NA),
                })
    return pd.DataFrame(rows, columns=["Source", "PathwayID", "Pathway", "Gene", "padj"])
