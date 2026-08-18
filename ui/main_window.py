from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QStackedWidget,
)

from ui.home_page import HomePage
from ui.page_16s import Page16S
from ui.page_feed_formula import PageFeedFormula
from ui.page_stats import PageStats
from ui.page_gas_fit import PageGasFit
from ui.page_silage_eval import PageSilageEval
from ui.page_rnaseq_analysis import PageRNASeqAnalysis
from ui.page_metabolomics import PageMetabolomics


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("田晨丰工具箱")
        self.resize(1280, 920)
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(8, 8, 8, 8)
        main_layout.setSpacing(8)

        nav_layout = QHBoxLayout()
        nav_layout.setSpacing(8)

        self.btn_home = QPushButton("首页")
        self.btn_page_16s = QPushButton("16S")
        self.btn_page_feed = QPushButton("饲料配方")
        self.btn_page_stats = QPushButton("统计分析")
        self.btn_page_gas_fit = QPushButton("体外产气拟合")
        self.btn_page_silage_eval = QPushButton("青贮综合评价")
        self.btn_page_rnaseq = QPushButton("RNASeq综合分析")

        self.btn_page_metabolomics = QPushButton("代谢组")

        self.btn_home.clicked.connect(self.show_home)
        self.btn_page_16s.clicked.connect(self.show_16s)
        self.btn_page_feed.clicked.connect(self.show_feed)
        self.btn_page_stats.clicked.connect(self.show_stats)
        self.btn_page_gas_fit.clicked.connect(self.show_gas_fit)
        self.btn_page_silage_eval.clicked.connect(self.show_silage_eval)
        self.btn_page_rnaseq.clicked.connect(self.show_rnaseq)

        self.btn_page_metabolomics.clicked.connect(self.show_metabolomics)

        nav_layout.addWidget(self.btn_home)
        nav_layout.addWidget(self.btn_page_16s)
        nav_layout.addWidget(self.btn_page_feed)
        nav_layout.addWidget(self.btn_page_stats)
        nav_layout.addWidget(self.btn_page_gas_fit)
        nav_layout.addWidget(self.btn_page_silage_eval)
        nav_layout.addWidget(self.btn_page_rnaseq)
        nav_layout.addWidget(self.btn_page_metabolomics)
        nav_layout.addStretch()

        self.stack = QStackedWidget()

        self.home_page = HomePage()
        self.page_16s = Page16S()
        self.page_feed_formula = PageFeedFormula()
        self.page_stats = PageStats()
        self.page_gas_fit = PageGasFit()
        self.page_silage_eval = PageSilageEval()
        self.page_rnaseq = PageRNASeqAnalysis()

        self.page_metabolomics = PageMetabolomics()

        self.home_page.on_open_16s = self.show_16s
        self.home_page.on_open_feed = self.show_feed
        self.home_page.on_open_stats = self.show_stats
        self.home_page.on_open_gas_fit = self.show_gas_fit
        self.home_page.on_open_silage_eval = self.show_silage_eval
        self.home_page.on_open_rnaseq = self.show_rnaseq

        self.home_page.on_open_metabolomics = self.show_metabolomics

        self.stack.addWidget(self.home_page)
        self.stack.addWidget(self.page_16s)
        self.stack.addWidget(self.page_feed_formula)
        self.stack.addWidget(self.page_stats)
        self.stack.addWidget(self.page_gas_fit)
        self.stack.addWidget(self.page_silage_eval)
        self.stack.addWidget(self.page_rnaseq)

        self.stack.addWidget(self.page_metabolomics)

        main_layout.addLayout(nav_layout)
        main_layout.addWidget(self.stack)

        self.setLayout(main_layout)

        self.apply_style()
        self.show_home()

    def show_home(self):
        self.stack.setCurrentIndex(0)

    def show_16s(self):
        self.stack.setCurrentIndex(1)

    def show_feed(self):
        self.stack.setCurrentIndex(2)

    def show_stats(self):
        self.stack.setCurrentIndex(3)

    def show_gas_fit(self):
        self.stack.setCurrentIndex(4)

    def show_silage_eval(self):
        self.stack.setCurrentIndex(5)

    def show_rnaseq(self):
        self.stack.setCurrentIndex(6)

    def show_metabolomics(self):
        self.stack.setCurrentWidget(self.page_metabolomics)

    def apply_style(self):
        self.setStyleSheet("""
            QWidget {
                background-color: #eef2f7;
                font-family: Microsoft YaHei, Arial;
                font-size: 11pt;
            }

            QPushButton {
                background-color: #4f46e5;
                color: white;
                border: none;
                border-radius: 8px;
                padding: 8px 14px;
                font-weight: bold;
                min-height: 24px;
            }

            QPushButton:hover {
                background-color: #4338ca;
            }

            QPushButton:pressed {
                background-color: #3730a3;
            }
        """)
