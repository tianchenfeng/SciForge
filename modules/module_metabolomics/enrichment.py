from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from modules.module_metabolomics.statistics import bh_fdr


def _split(value) -> list[str]:
    if pd.isna(value) or str(value).strip() in {"", "-", "nan"}:
        return []
    return [item.strip() for item in str(value).replace("|", ";").split(";") if item.strip() and item.strip() != "-"]


def build_pathway_mapping(features: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in features.iterrows():
        compounds = _split(row.get("KEGG", ""))
        pathways = _split(row.get("Pathway IDs", ""))
        names = _split(row.get("Pathway names", ""))
        for index, pathway in enumerate(pathways):
            pathway_name = names[index] if index < len(names) else pathway
            for compound in compounds:
                rows.append({"Pathway ID": pathway, "Pathway": pathway_name, "KEGG": compound, "Feature ID": str(row.get("Feature ID", "")), "Metabolite": str(row.get("Metabolite name", ""))})
    return pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame(columns=["Pathway ID", "Pathway", "KEGG", "Feature ID", "Metabolite"])


def run_pathway_enrichment(significant_features, features: pd.DataFrame, mapping: pd.DataFrame | None = None) -> pd.DataFrame:
    mapping = build_pathway_mapping(features) if mapping is None else mapping.copy()
    if mapping.empty:
        return pd.DataFrame(columns=["Pathway", "Pathway ID", "Hit count", "Background count", "Study compounds", "Background compounds", "Enrichment ratio", "P", "FDR", "Hit metabolites"])
    significant = set(map(str, significant_features))
    study = set(mapping.loc[mapping["Feature ID"].astype(str).isin(significant), "KEGG"].astype(str))
    universe = set(mapping["KEGG"].astype(str))
    rows = []
    for (pathway_id, pathway), block in mapping.groupby(["Pathway ID", "Pathway"], dropna=False):
        background = set(block["KEGG"].astype(str)); hits = background & study
        if not hits: continue
        p = hypergeom.sf(len(hits) - 1, len(universe), len(background), len(study))
        names = block[block["KEGG"].astype(str).isin(hits)]["Metabolite"].replace("", pd.NA).dropna().astype(str).unique()
        rows.append({"Pathway": pathway, "Pathway ID": pathway_id, "Hit count": len(hits), "Background count": len(background), "Study compounds": len(study), "Background compounds": len(universe), "Enrichment ratio": len(hits) / max(len(study), 1), "P": p, "Hit metabolites": "; ".join(names)})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["FDR"] = bh_fdr(result["P"])
        result = result.sort_values(["FDR", "P"]).reset_index(drop=True)
    return result

