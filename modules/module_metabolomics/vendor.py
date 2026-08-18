from __future__ import annotations

from pathlib import Path

import pandas as pd

from modules.module_metabolomics.adapters.base import read_table_smart


def load_vendor_oplsda(dataset) -> dict[str, pd.DataFrame]:
    root_text = dataset.vendor_results.get("oplsda_root")
    if not root_text:
        return {}
    root = Path(root_text); outputs = {}
    for comparison_dir in sorted(path for path in root.iterdir() if path.is_dir()):
        for key, name in {"model": "OPLS-DA.model.xls", "scores": "OPLS-DA.sites.xls", "vip": "OPLS-DA.vips.xls", "permutation": "OPLS-DA.permMN.xls"}.items():
            path = comparison_dir / name
            if path.exists():
                table = read_table_smart(path, low_memory=False)
                table.insert(0, "Comparison", comparison_dir.name)
                table.insert(1, "Source", "Upstream vendor result; not recalculated by SciForge")
                outputs.setdefault(key, []).append(table)
    return {key: pd.concat(tables, ignore_index=True) for key, tables in outputs.items()}

