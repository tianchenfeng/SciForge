import json
import os

import numpy as np
import pandas as pd
from PySide6.QtWidgets import QApplication

from modules.module_metabolomics.service import AnalysisConfig, MetabolomicsAnalysisService
from ui.page_metabolomics import PageMetabolomics


def test_end_to_end_result_package(tmp_path):
    rng = np.random.default_rng(31); samples = [f"S{i}" for i in range(12)] + ["QC1", "QC2", "QC3"]
    matrix = rng.lognormal(4, .2, (25, 15)); matrix[:3, 6:12] *= 2
    frame = pd.DataFrame(matrix, columns=samples); frame.insert(0, "Feature ID", [f"F{i}" for i in range(25)])
    input_path = tmp_path / "matrix.csv"; frame.to_csv(input_path, index=False)
    records = []
    for sample in samples:
        is_qc = sample.startswith("QC")
        records.append({"Sample ID": sample, "Sample type": "qc" if is_qc else "biological", "Group": "QC" if is_qc else ("A" if int(sample[1:]) < 6 else "B"), "Factor A": "", "Factor B": ""})
    config = AnalysisConfig(input_path=str(input_path), output_parent=str(tmp_path), metadata_records=records, confirmed_input_scale="raw", permutations=3, heatmap_top_n=10)
    result = MetabolomicsAnalysisService.run(config)
    assert result["success"] is True
    output = tmp_path / next(path.name for path in tmp_path.iterdir() if path.is_dir() and path.name.startswith("Metabolomics_Analysis_"))
    assert (output / "tables" / "metabolomics_analysis.xlsx").exists()
    assert (output / "report" / "metabolomics_report.html").exists()
    manifest = json.loads((output / "00_manifest" / "analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["random seed"] == 2026
    assert manifest["fold_change"]["fc_formula"] == "mean(target) / mean(reference)"
    assert manifest["PLSDA"]["Scaling fitted within folds"] is True
    assert "not a universal" in manifest["PLSDA"]["Q2 interpretation rules"]["scope"]
    assert "Significant features" not in result["summary"]
    assert "Significant feature-contrast records" in result["summary"]
    summary_text = (output / "report" / "analysis_summary.md").read_text(encoding="utf-8")
    html_text = (output / "report" / "metabolomics_report.html").read_text(encoding="utf-8")
    interpretation = result["summary"]["PLS-DA interpretation"]
    assert interpretation in summary_text and interpretation in html_text
    excel_summary = pd.read_excel(output / "tables" / "metabolomics_analysis.xlsx", sheet_name="Summary")
    excel_values = dict(zip(excel_summary["Metric"], excel_summary["Value"]))
    assert excel_values["PLS-DA interpretation"] == interpretation


def test_gui_smoke(monkeypatch):
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication([])
    page = PageMetabolomics()
    assert page.results_tabs.count() >= 13
    assert page.scale_combo.currentData() == "unknown"
    page.close(); app.processEvents()
