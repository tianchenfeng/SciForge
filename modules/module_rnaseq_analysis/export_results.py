from __future__ import annotations

from pathlib import Path

import pandas as pd


REQUIRED_SHEETS = [
    "all_deg",
    "filtered_deg",
    "up_genes",
    "down_genes",
    "deg_with_annotation",
    "filtered_deg_with_annotation",
    "up_genes_with_annotation",
    "down_genes_with_annotation",
    "candidate_genes_ranked",
    "top_GO",
    "significant_GO",
    "keyword_GO",
    "top_KEGG",
    "significant_KEGG",
    "keyword_KEGG",
    "pathway_gene_map",
]


def export_excel(tables: dict[str, pd.DataFrame], output_path: str | Path) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for sheet_name in REQUIRED_SHEETS:
            tables.get(sheet_name, pd.DataFrame()).to_excel(writer, index=False, sheet_name=sheet_name)
    return str(output_path)
