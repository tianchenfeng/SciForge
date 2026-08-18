from __future__ import annotations

import base64
import hashlib
import html
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SHEET_ORDER = ["Summary", "Samples", "QC_Summary", "QC", "Feature_Metadata", "Annotations", "PCA", "PCA_Loadings", "PCA_Variance", "PLSDA", "PLSDA_OOF", "PLSDA_VIP", "PLSDA_CV", "PLSDA_Permutation", "Overall_Test", "Contrasts", "TwoFactor_A", "TwoFactor_B", "TwoFactor_Interaction", "SimpleEffects", "Differential", "Intersections", "Pathways", "Correlations", "Parameters", "Vendor_OPLS_Model", "Vendor_OPLS_Scores", "Vendor_OPLS_VIP", "Vendor_OPLS_Permutation"]


def file_sha256(path: str | Path) -> str | None:
    path = Path(path)
    if not path.is_file(): return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict): return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [json_safe(v) for v in value]
    if isinstance(value, Path): return str(value)
    if isinstance(value, (np.integer,)): return int(value)
    if isinstance(value, (np.floating,)): return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.ndarray): return json_safe(value.tolist())
    if pd.isna(value) if not isinstance(value, (str, bool)) else False: return None
    return value


def export_excel(tables: dict[str, pd.DataFrame], path: str | Path) -> str:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        written = False
        for name in SHEET_ORDER:
            table = tables.get(name)
            if table is None: continue
            if not isinstance(table, pd.DataFrame): table = pd.DataFrame(table)
            table.to_excel(writer, sheet_name=name[:31], index=True if name in {"QC"} and table.index.name else False)
            written = True
        if not written: pd.DataFrame({"Message": ["No result tables were produced."]}).to_excel(writer, sheet_name="Summary", index=False)
    return str(path)


def export_manifest(manifest: dict, path: str | Path) -> str:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(manifest), ensure_ascii=False, indent=2), encoding="utf-8")
    return str(path)


def _image_data(path: Path) -> str:
    mime = "image/svg+xml" if path.suffix.lower() == ".svg" else "image/png"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"


def export_html_report(summary: dict, sections: list[tuple[str, str]], figure_paths: list[str], warnings: list[str], path: str | Path) -> str:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    summary_rows = "".join(f"<tr><th>{html.escape(str(k))}</th><td>{html.escape(str(v))}</td></tr>" for k, v in summary.items())
    section_html = "".join(f"<section><h2>{html.escape(title)}</h2><p>{html.escape(text)}</p></section>" for title, text in sections)
    images = []
    for item in figure_paths:
        image = Path(item)
        if image.exists() and image.suffix.lower() == ".png": images.append(f'<figure><img src="{_image_data(image)}"><figcaption>{html.escape(image.stem)}</figcaption></figure>')
    warning_html = "".join(f"<li>{html.escape(str(w))}</li>" for w in warnings) or "<li>None</li>"
    document = f'''<!doctype html><html><head><meta charset="utf-8"><title>SciForge Metabolomics Report</title><style>body{{font-family:Arial,'Microsoft YaHei',sans-serif;max-width:1100px;margin:auto;padding:30px;color:#1f2937}}h1,h2{{color:#4338ca}}table{{border-collapse:collapse;width:100%}}th,td{{border-bottom:1px solid #ddd;padding:8px;text-align:left}}figure{{margin:28px 0}}img{{max-width:100%;height:auto}}.note{{background:#fef3c7;padding:12px;border-left:4px solid #d97706}}</style></head><body><h1>SciForge Metabolomics Analysis</h1><p class="note">Statistical facts are reported separately from biological or causal interpretation. VIP is not a significance test.</p><table>{summary_rows}</table>{section_html}<h2>Warnings</h2><ul>{warning_html}</ul><h2>Figures</h2>{''.join(images)}</body></html>'''
    path.write_text(document, encoding="utf-8")
    return str(path)


def write_text(path: str | Path, content: str) -> str:
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True); path.write_text(content, encoding="utf-8"); return str(path)
