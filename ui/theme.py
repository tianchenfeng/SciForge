"""Shared desktop appearance and layout safeguards, independent of analysis."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox, QLabel, QLayout, QPushButton, QSizePolicy, QTableWidget,
    QTabWidget, QTextEdit, QWidget,
)

APP_STYLE = """
QWidget { font-family: 'Microsoft YaHei UI', 'Segoe UI'; font-size: 13px; color: #24272d; }
QWidget#appWindow { background: #f5f5f7; }
QFrame#sidebar { background: #eceef2; border-right: 1px solid #dedfe5; }
QLabel { background: transparent; }
QLabel#brandMark { color: white; font-size: 28px; font-weight: 700;
    border-radius: 12px; background: qlineargradient(x1:0,y1:0,x2:1,y2:1,stop:0 #65aaff,stop:1 #3664e8); }
QLabel#brandName { font-size: 25px; font-weight: 700; color: #1d1d1f; margin-top: 5px; }
QLabel#mutedLabel, QLabel#subtitleLabel, QLabel#hintLabel { color: #737783; }
QLabel#navCaption { color: #858994; font-size: 11px; margin: 0 9px 6px 9px; }
QLabel#sidebarFooter { color: #8b8e98; font-size: 10px; }
QLabel#workspaceStatus { color: #598171; font-size: 11px; }
QLabel#pageTitle { font-size: 24px; font-weight: 700; color: #1d1d1f; margin-bottom: 10px; }
QPushButton { background: white; color: #404752; border: 1px solid #dce0e6;
    border-radius: 8px; padding: 8px 13px; min-height: 18px; }
QPushButton:hover { background: #eef4ff; border-color: #a6c7ff; color: #0068d9; }
QPushButton:pressed { background: #dfeaff; }
QPushButton:focus { border: 1px solid #007aff; }
QPushButton:disabled { background: #eff0f3; color: #a8abb4; border-color: #e4e6eb; }
QPushButton[primary="true"] { background: #007aff; color: white; border-color: #007aff; font-weight: 600; }
QPushButton[primary="true"]:hover { background: #0067df; }
QPushButton[primary="true"]:disabled { background: #b6d5fa; color: #f8fbff; border-color: #b6d5fa; }
QPushButton#navButton { text-align: left; background: transparent; border: 1px solid transparent;
    border-radius: 9px; padding: 10px 14px; min-height: 18px; }
QPushButton#navButton:hover { background: #e1e5ed; color: #24272d; }
QPushButton#navButton:checked { background: #dce8fc; color: #0066d5; font-weight: 600; }
QPushButton[compact="true"] { padding: 0; min-height: 0; border-radius: 7px; font-size: 18px; }
QGroupBox { background: white; border: 1px solid #e1e4ea; border-radius: 12px;
    margin-top: 20px; padding: 18px 12px 12px 12px; font-weight: 600; }
QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left;
    left: 14px; padding: 0 5px; color: #5c6370; }
QGroupBox QGroupBox { background: #fafbfc; }
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox { background: #f9fafc;
    border: 1px solid #dfe3e9; border-radius: 7px; padding: 7px 9px; min-height: 19px; selection-background-color: #c8dfff; }
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus { border-color: #007aff; background: white; }
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled { color: #a1a6b0; background: #f0f1f4; }
QComboBox { padding-right: 28px; }
QComboBox::drop-down { width: 23px; border: none; }
QComboBox QAbstractItemView { background: white; selection-background-color: #dceaff; selection-color: #0066d5; }
QTextEdit, QTextBrowser, QTableWidget { background: white; border: 1px solid #e1e4ea;
    border-radius: 9px; padding: 6px; selection-background-color: #dceaff; selection-color: #1d1d1f; }
QTableWidget { gridline-color: #edf0f4; alternate-background-color: #f8faff; }
QHeaderView::section { background: #f3f5f8; border: none; border-bottom: 1px solid #e1e4ea; padding: 9px 7px; color: #636b78; }
QCheckBox, QRadioButton { spacing: 8px; padding: 4px 0; background: transparent; }
QCheckBox::indicator, QRadioButton::indicator { width: 16px; height: 16px; }
QProgressBar { background: #e8ecf2; border: none; border-radius: 5px; min-height: 10px; max-height: 18px; text-align: center; font-size: 10px; }
QProgressBar::chunk { background: #007aff; border-radius: 5px; }
QScrollArea#pageScroll, QScrollArea#controlsScroll, QScrollArea > QWidget > QWidget { background: transparent; border: none; }
QScrollBar:vertical { background: transparent; width: 9px; margin: 2px; }
QScrollBar::handle:vertical { background: #c7cbd3; border-radius: 3px; min-height: 32px; }
QScrollBar:horizontal { background: transparent; height: 9px; margin: 2px; }
QScrollBar::handle:horizontal { background: #c7cbd3; border-radius: 3px; min-width: 32px; }
QScrollBar::add-line, QScrollBar::sub-line { width: 0; height: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }
QTabWidget::pane { background: white; border: 1px solid #e1e4ea; border-radius: 9px; }
QTabBar::tab { background: transparent; color: #7c8290; padding: 9px 12px; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { color: #007aff; border-bottom-color: #007aff; }
QSplitter::handle { background: #e5e7ec; width: 5px; height: 5px; }
QFrame#hero { border: 1px solid #e1e7f1; border-radius: 20px;
    background: qlineargradient(x1:0,y1:0,x2:1,y2:1,stop:0 #ffffff,stop:.55 #f1f6ff,stop:1 #e8edfc); }
QLabel#eyebrow { color: #006cdb; font-size: 10px; font-weight: 600; }
QLabel#heroTitle { color: #1d1d1f; font-size: 34px; font-weight: 700; }
QLabel#heroSubtitle { color: #737c8d; font-size: 13px; }
QLabel#sectionTitle { font-size: 18px; font-weight: 600; }
QPushButton#moduleCard { background: white; border: 1px solid #e2e5eb; border-radius: 14px; padding: 0; min-height: 142px; text-align: left; }
QPushButton#moduleCard:hover { background: #f9fbff; border-color: #9cbff1; }
QLabel#cardTitle { font-size: 15px; font-weight: 600; color: #303744; }
QLabel#cardDescription { font-size: 11px; color: #8b919e; }
QLabel#cardArrow { color: #a3adbd; font-size: 18px; }
"""


def prepare_page(page):
    """Remove conflicting legacy styles and allow long content to scroll safely."""
    for widget in [page, *page.findChildren(QWidget)]:
        widget.setStyleSheet("")
        if isinstance(widget, QLabel):
            widget.setWordWrap(True)
            widget.setMinimumWidth(0)
            widget.setTextFormat(Qt.PlainText)
        if isinstance(widget, QComboBox):
            widget.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
            widget.setMinimumContentsLength(8)
        if isinstance(widget, QTextEdit):
            widget.setMinimumHeight(140)
        if isinstance(widget, QTableWidget):
            widget.setMinimumHeight(200)
            widget.setAlternatingRowColors(True)
            widget.verticalHeader().setDefaultSectionSize(40)
        if isinstance(widget, QTabWidget):
            widget.setUsesScrollButtons(True)
            widget.tabBar().setExpanding(False)
        if isinstance(widget, QPushButton):
            widget.setCursor(Qt.PointingHandCursor)
            if widget.text() in {"开始分析", "开始计算", "识别并校验数据"}:
                widget.setProperty("primary", True)
    first = page.layout().itemAt(0).widget()
    if isinstance(first, QLabel):
        first.setObjectName("pageTitle")
        first.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
    page.layout().setSizeConstraint(QLayout.SetMinimumSize)
    page.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
