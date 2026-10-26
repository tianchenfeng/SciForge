"""Optional, non-interactive diagnostics for a distributed Windows build."""
import argparse
import importlib
import json
import os
from pathlib import Path
import sys
import traceback


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = []

    def check(name, fn):
        try:
            detail = fn()
            checks.append({"name":name,"passed":True,"detail":detail})
        except Exception:
            checks.append({"name":name,"passed":False,"error":traceback.format_exc()})

    for name in ["PySide6.QtWidgets", "numpy", "pandas", "openpyxl", "xlrd", "docx", "scipy.stats", "statsmodels.api",
                 "statsmodels.formula.api", "sklearn.decomposition", "sklearn.cross_decomposition", "seaborn", "networkx",
                 "matplotlib.backends.backend_agg", "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg"]:
        check("import:"+name, lambda n=name: str(importlib.import_module(n).__name__))

    def ui_check():
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from PySide6.QtWidgets import QApplication
        from PySide6.QtGui import QFontDatabase
        from ui.main_window import MainWindow
        app = QApplication.instance() or QApplication([])
        if not QFontDatabase.families():
            for name in ["msyh.ttc", "msyhbd.ttc", "segoeui.ttf"]:
                QFontDatabase.addApplicationFont("C:/Windows/Fonts/"+name)
        window = MainWindow()
        window.show()
        snapshots = []
        for width, height in [(800,560),(1360,940)]:
            window.resize(width,height)
            for index,page in enumerate(window.pages):
                window._show_page(index)
                for _ in range(8): app.processEvents()
                assert window.size().toTuple() == (width,height)
                assert window.page_scrolls[index].horizontalScrollBar().maximum() == 0
                snapshots.append([index,width,height])
                if index in (0,3): window.grab().save(str(output/f"frozen-page-{index}-{width}.png"))
        window.close()
        return snapshots
    check("eight-page-layout", ui_check)

    def stats_check():
        import numpy as np
        import pandas as pd
        from modules.module_stats.service import run_one_way_batch_analysis, run_two_way_batch_analysis
        frame = pd.DataFrame([{"Group":g,"Time":t,"value":10+gi*2+ti+r*.2}
                              for gi,g in enumerate(["C","T"]) for ti,t in enumerate(["D1","D2"]) for r in range(3)])
        for mode,service in [("one",run_one_way_batch_analysis),("two",run_two_way_batch_analysis)]:
            source = output/f"{mode}.csv"
            data = frame.drop(columns="Time") if mode=="one" else frame.copy()
            data.to_csv(source,index=False)
            result = service(source,output/mode)
            assert result["success"],result["message"]
            assert Path(result["word_path"]).is_file() and Path(result["excel_path"]).is_file()
            data.loc[0,"value"] = np.inf
            data.to_csv(source,index=False)
            assert service(source,output/mode)["success"] is False
        return "one/two-way exports and infinity rejection passed"
    check("statistics-and-exports", stats_check)
    base = Path(getattr(sys,"_MEIPASS",Path(__file__).resolve().parents[1]))
    check("embedded-16s-scripts", lambda: _check_scripts(base))
    result = {"success":all(c["passed"] for c in checks),"frozen":bool(getattr(sys,"frozen",False)),"checks":checks}
    (output/"runtime-check.json").write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding="utf-8")
    return 0 if result["success"] else 1


def _check_scripts(base):
    from modules.module_16s.config import PIPELINE_ORDER
    for name in PIPELINE_ORDER:
        assert (base/"modules/module_16s/scripts"/name).is_file(),name
    return PIPELINE_ORDER
