from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication, QButtonGroup, QFrame, QHBoxLayout, QLabel, QPushButton,
    QScrollArea, QStackedWidget, QVBoxLayout, QWidget,
)

from ui.home_page import HomePage
from ui.page_16s import Page16S
from ui.page_feed_formula import PageFeedFormula
from ui.page_stats import PageStats
from ui.page_gas_fit import PageGasFit
from ui.page_silage_eval import PageSilageEval
from ui.page_rnaseq_analysis import PageRNASeqAnalysis
from ui.page_metabolomics import PageMetabolomics
from ui.theme import APP_STYLE, prepare_page


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setObjectName("appWindow")
        self.setWindowTitle("SciForge · 科研数据工作台")
        self.setMinimumSize(800, 560)
        screen = QApplication.primaryScreen().availableGeometry()
        self.resize(min(1360, screen.width()), min(940, screen.height()))
        self.init_ui()

    def init_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        self.sidebar = QFrame()
        self.sidebar.setObjectName("sidebar")
        self.sidebar.setFixedWidth(206)
        nav = QVBoxLayout(self.sidebar)
        nav.setContentsMargins(16, 26, 16, 20)
        nav.setSpacing(7)
        brand = QLabel("S")
        brand.setObjectName("brandMark")
        brand.setFixedSize(42, 42)
        brand.setAlignment(Qt.AlignCenter)
        nav.addWidget(brand)
        name = QLabel("SciForge")
        name.setObjectName("brandName")
        nav.addWidget(name)
        subtitle = QLabel("科研数据工作台")
        subtitle.setObjectName("mutedLabel")
        nav.addWidget(subtitle)
        nav.addSpacing(28)
        caption = QLabel("工作空间")
        caption.setObjectName("navCaption")
        nav.addWidget(caption)
        self.nav_group = QButtonGroup(self)
        self.nav_group.setExclusive(True)
        names = ["总览", "16S 微生物组", "饲料配方", "统计分析", "体外产气拟合",
                 "青贮综合评价", "RNA-Seq 分析", "代谢组分析"]
        attributes = ["btn_home", "btn_page_16s", "btn_page_feed", "btn_page_stats",
                      "btn_page_gas_fit", "btn_page_silage_eval", "btn_page_rnaseq",
                      "btn_page_metabolomics"]
        for index, (text, attribute) in enumerate(zip(names, attributes)):
            button = QPushButton(text)
            button.setObjectName("navButton")
            button.setCheckable(True)
            button.setCursor(Qt.PointingHandCursor)
            button.clicked.connect(lambda checked=False, i=index: self._show_page(i))
            self.nav_group.addButton(button, index)
            setattr(self, attribute, button)
            nav.addWidget(button)
        nav.addStretch()
        footer = QLabel("为每一次探索而造\nCreated by Tian Chenfeng")
        footer.setWordWrap(True)
        footer.setObjectName("sidebarFooter")
        nav.addWidget(footer)
        root.addWidget(self.sidebar)
        body = QWidget()
        body_layout = QVBoxLayout(body)
        body_layout.setContentsMargins(22, 18, 22, 16)
        body_layout.setSpacing(12)
        header = QHBoxLayout()
        self.breadcrumb = QLabel("工作空间 / 总览")
        self.breadcrumb.setObjectName("mutedLabel")
        header.addWidget(self.breadcrumb)
        header.addStretch()
        status = QLabel("●  本地工作空间")
        status.setObjectName("workspaceStatus")
        header.addWidget(status)
        body_layout.addLayout(header)
        self.stack = QStackedWidget()
        self.stack.setMinimumSize(0, 0)
        body_layout.addWidget(self.stack, 1)
        root.addWidget(body, 1)
        self.home_page = HomePage()
        self.page_16s = Page16S()
        self.page_feed_formula = PageFeedFormula()
        self.page_stats = PageStats()
        self.page_gas_fit = PageGasFit()
        self.page_silage_eval = PageSilageEval()
        self.page_rnaseq = PageRNASeqAnalysis()
        self.page_metabolomics = PageMetabolomics()
        self.pages = [self.home_page, self.page_16s, self.page_feed_formula, self.page_stats,
                      self.page_gas_fit, self.page_silage_eval, self.page_rnaseq, self.page_metabolomics]
        self.page_scrolls = []
        for index, page in enumerate(self.pages):
            if index:
                prepare_page(page)
            scroll = QScrollArea()
            scroll.setObjectName("pageScroll")
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.NoFrame)
            scroll.setWidget(page)
            scroll.setMinimumSize(0, 0)
            self.page_scrolls.append(scroll)
            self.stack.addWidget(scroll)
        for callback, handler in [
            ("on_open_16s", self.show_16s), ("on_open_feed", self.show_feed),
            ("on_open_stats", self.show_stats), ("on_open_gas_fit", self.show_gas_fit),
            ("on_open_silage_eval", self.show_silage_eval), ("on_open_rnaseq", self.show_rnaseq),
            ("on_open_metabolomics", self.show_metabolomics),
        ]:
            setattr(self.home_page, callback, handler)
        self.page_names = names
        self.apply_style()
        self.show_home()

    def _show_page(self, index):
        self.stack.setCurrentIndex(index)
        self.nav_group.button(index).setChecked(True)
        self.breadcrumb.setText(f"工作空间 / {self.page_names[index]}")

    def show_home(self): self._show_page(0)
    def show_16s(self): self._show_page(1)
    def show_feed(self): self._show_page(2)
    def show_stats(self): self._show_page(3)
    def show_gas_fit(self): self._show_page(4)
    def show_silage_eval(self): self._show_page(5)
    def show_rnaseq(self): self._show_page(6)
    def show_metabolomics(self): self._show_page(7)

    def apply_style(self):
        self.setStyleSheet(APP_STYLE)
