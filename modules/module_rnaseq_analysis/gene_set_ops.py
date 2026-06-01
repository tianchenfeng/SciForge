from __future__ import annotations

import re

import pandas as pd


TEXT_COLUMNS = ["description", "Description", "pathway", "Pathway", "GO", "KEGG", "Term", "ID"]


def parse_keywords(keyword_text: str | None) -> list[str]:
    if not keyword_text:
        return []
    normalized = keyword_text.replace(";", ",").replace("\n", ",")
    return [item.strip() for item in normalized.split(",") if item.strip()]


def keyword_hit_count(row: pd.Series, keywords: list[str]) -> int:
    if not keywords:
        return 0
    text = " ".join(str(row[col]) for col in TEXT_COLUMNS if col in row.index and pd.notna(row[col])).lower()
    return sum(1 for keyword in keywords if keyword.lower() in text)


def split_gene_list(value) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in re.split(r"[/,;|]\s*", str(value)) if item.strip()]
