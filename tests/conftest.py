"""Desktop acceptance tests keep every generated input/output in the local test folder."""
import os
from pathlib import Path
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pytest


@pytest.fixture(scope="session")
def app():
    from PySide6.QtWidgets import QApplication
    from PySide6.QtGui import QFontDatabase
    instance = QApplication.instance() or QApplication([])
    # The Windows offscreen plugin does not enumerate installed system fonts.
    if not QFontDatabase.families():
        for name in ("msyh.ttc", "msyhbd.ttc", "segoeui.ttf"):
            font = Path("C:/Windows/Fonts") / name
            if font.exists():
                QFontDatabase.addApplicationFont(str(font))
    return instance
