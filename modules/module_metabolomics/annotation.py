from __future__ import annotations

import re

import numpy as np
import pandas as pd


CAUTION_RULES = {
    "Antibiotic-like annotation": r"cephalospor|penicillin|antibiotic|tetracyclin|macrolide",
    "Drug-like annotation": r"drug|pharmaceutical|statin|ibuprofen|acetaminophen",
    "Xenobiotic-like annotation": r"xenobiotic|pesticide|herbicide|synthetic",
    "Plant hormone-like annotation": r"gibberellin|auxin|cytokinin|abscisic acid|brassinosteroid",
}


def audit_annotations(features: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    result = features.copy()
    for column in ["Feature ID", "Ion mode", "m/z", "RT", "Metabolite name", "Formula", "KEGG", "HMDB", "Library ID", "Annotation level", "MS/MS score"]:
        if column not in result:
            result[column] = ""
    level = result["Annotation level"].fillna("").astype(str)
    identified = ~level.isin(["", "-", "unknown", "nan"])
    result["Identified"] = identified
    result["Entity type"] = np.where(identified, "Annotated metabolite", "Feature")
    text = result["Metabolite name"].fillna("").astype(str)
    cautions = []
    for index, (value, is_identified) in enumerate(zip(text, identified)):
        flags = [label for label, pattern in CAUTION_RULES.items() if is_identified and re.search(pattern, value, flags=re.I)]
        if is_identified and str(value).strip() and not flags and str(level.iloc[index]) in {"B(ii)", "Bii", "Level 3", "Level 4"}:
            flags.append("Low-confidence annotation")
        cautions.append("; ".join(flags))
    result["Annotation caution"] = cautions
    coverage = pd.DataFrame({"Status": ["Annotated", "Unannotated"], "Count": [int(identified.sum()), int((~identified).sum())]})
    confidence = level.where(identified, "Unannotated").value_counts().rename_axis("Annotation level").reset_index(name="Count")
    ion_mode = result["Ion mode"].fillna("unknown").replace("", "unknown").value_counts().rename_axis("Ion mode").reset_index(name="Count")
    return result, {"coverage": coverage, "confidence": confidence, "ion_mode": ion_mode}
