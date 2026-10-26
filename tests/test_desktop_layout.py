from itertools import combinations
import json
from pathlib import Path

import pytest
from PySide6.QtCore import QPoint, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QAbstractButton, QCheckBox, QComboBox, QLineEdit, QPushButton, QRadioButton, QScrollArea, QTabWidget, QWidget

from ui.main_window import MainWindow

OUT = Path(__file__).resolve().parents[1] / "test" / "acceptance_20261026"


@pytest.fixture
def window(app):
    widget = MainWindow()
    widget.show()
    yield widget
    widget.close()
    widget.deleteLater()
    app.processEvents()


def settle(app):
    for _ in range(12):
        app.processEvents()


@pytest.mark.parametrize("size", [(800, 560), (960, 640), (1024, 768), (1280, 720), (1360, 940), (1920, 1080)])
def test_every_page_geometry(app, window, size):
    OUT.mkdir(parents=True, exist_ok=True)
    window.resize(*size)
    settle(app)
    report = []
    for index, page in enumerate(window.pages):
        window.nav_group.button(index).click()
        settle(app)
        assert window.stack.currentIndex() == index
        assert window.nav_group.button(index).isChecked()
        assert window.size().toTuple() == size
        # All actionable controls must be inside their immediate content parent.
        controls = [w for w in page.findChildren(QWidget)
                    if isinstance(w, (QPushButton, QCheckBox, QRadioButton, QLineEdit, QComboBox)) and w.isVisibleTo(page)]
        for control in controls:
            parent = control.parentWidget()
            if parent.objectName() == "qt_scrollarea_viewport":
                continue  # Table cells are clipped and scrolled by the table viewport.
            rect = control.geometry()
            assert rect.width() > 0 and rect.height() > 0, (index, type(control).__name__, control.objectName(), parent, rect)
            assert parent.rect().adjusted(-1, -1, 1, 1).contains(rect), (
                index, size, control.objectName(), getattr(control, "text", lambda: "")(), rect, parent.rect())
        # Sibling controls cannot cover one another; Qt's own spinbox editor is excluded.
        for first, second in combinations(controls, 2):
            if first.parentWidget() is second.parentWidget():
                assert not first.geometry().intersects(second.geometry()), (index, size, first, second)
        for button in window.nav_group.buttons():
            assert window.sidebar.rect().contains(button.geometry())
        for first, second in combinations(window.nav_group.buttons(), 2):
            assert not first.geometry().intersects(second.geometry())
        scroll = window.page_scrolls[index]
        # Every button is reachable by scrolling, even when it is below the fold.
        for button in [w for w in controls if isinstance(w, QAbstractButton)]:
            ancestor = button.parentWidget()
            while ancestor is not window:
                if isinstance(ancestor, QScrollArea):
                    ancestor.ensureWidgetVisible(button, 8, 8)
                ancestor = ancestor.parentWidget()
            settle(app)
            assert window.rect().contains(button.mapTo(window, button.rect().center())), (index, size, button)
        scroll.verticalScrollBar().setValue(0)
        scroll.horizontalScrollBar().setValue(0)
        if size in [(800, 560), (1360, 940)]:
            window.grab().save(str(OUT / f"page-{index}-{size[0]}x{size[1]}.png"))
        report.append({"page": index, "content": page.size().toTuple(), "controls": len(controls),
                       "horizontal_scroll": scroll.horizontalScrollBar().maximum(),
                       "vertical_scroll": scroll.verticalScrollBar().maximum()})
    (OUT / f"layout-{size[0]}x{size[1]}.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


def test_home_shortcuts_and_mode_switch(app, window):
    destinations = {"stats": 3, "16s": 1, "metabolomics": 7, "rnaseq": 6,
                    "gas_fit": 4, "silage_eval": 5, "feed": 2}
    for key, index in destinations.items():
        window.show_home()
        getattr(window.home_page, f"btn_{key}").click()
        assert window.stack.currentIndex() == index
    window.show_stats()
    stats = window.page_stats
    stats.combo_mode.setCurrentIndex(1)
    assert not stats.combo_factor_b.isHidden()
    stats.combo_mode.setCurrentIndex(0)
    assert stats.combo_factor_b.isHidden()
    assert len(window.page_metabolomics.results_tabs.findChildren(QTabWidget)) == 0


def test_feed_constraint_buttons_act_on_the_displayed_values(app, window):
    window.show_feed()
    settle(app)
    for widgets in window.page_feed_formula.constraint_widgets.values():
        widgets["enable"].setChecked(True)
        for key in ("target", "error"):
            spin = widgets[key]
            before = spin.value()
            QTest.mouseClick(widgets[f"{key}_plus"], Qt.LeftButton)
            assert spin.value() == pytest.approx(before + spin.singleStep())
            QTest.mouseClick(widgets[f"{key}_minus"], Qt.LeftButton)
            assert spin.value() == pytest.approx(before)
            spin.setValue(spin.maximum())
            widgets[f"{key}_plus"].click()
            assert spin.value() == spin.maximum()
            spin.setValue(0)
            widgets[f"{key}_minus"].click()
            assert spin.value() == 0
        widgets["enable"].setChecked(False)
        for key in ("target", "error", "target_plus", "target_minus", "error_plus", "error_minus"):
            assert not widgets[key].isEnabled()


def test_long_paths_and_column_names(app, window):
    path = "F:/" + "毕业实验数据_很长的文件夹/" * 18 + "单因素双因素数据.xlsx"
    window.page_16s.input_label.setText(path)
    window.page_16s.output_label.setText(path)
    window.page_stats.set_columns(["因素_" + "复杂名称" * 30, "时间", "指标"])
    window.page_stats.input_edit.setText(path)
    window.page_metabolomics.input_edit.setText(path)
    window.resize(800, 560)
    for index in (1, 3, 7):
        window._show_page(index)
        settle(app)
        assert window.width() == 800
        assert window.page_scrolls[index].viewport().width() > 400
