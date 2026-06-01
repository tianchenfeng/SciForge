from PySide6.QtWidgets import QWidget, QVBoxLayout, QLabel, QPushButton
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont


class HomePage(QWidget):
    def __init__(self):
        super().__init__()
        self.on_open_16s = None
        self.on_open_feed = None
        self.on_open_stats = None
        self.on_open_gas_fit = None
        self.on_open_silage_eval = None
        self.on_open_rnaseq = None
        self.init_ui()
        self.apply_style()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(40, 40, 40, 40)
        layout.setSpacing(24)

        title = QLabel("田晨丰工具箱")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Microsoft YaHei", 20, QFont.Bold))

        subtitle = QLabel("请选择要使用的功能模块")
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setObjectName("subtitleLabel")

        layout.addWidget(title)
        layout.addWidget(subtitle)

        btn_layout = QVBoxLayout()
        btn_layout.setSpacing(18)

        self.btn_16s = QPushButton("16S 生信分析")
        self.btn_feed = QPushButton("智能饲料配方")
        self.btn_stats = QPushButton("通用统计分析")
        self.btn_gas_fit = QPushButton("体外产气拟合")
        self.btn_silage_eval = QPushButton("青贮综合评价")
        self.btn_rnaseq = QPushButton("RNASeq 综合分析")

        self.btn_16s.clicked.connect(self.handle_open_16s)
        self.btn_feed.clicked.connect(self.handle_open_feed)
        self.btn_stats.clicked.connect(self.handle_open_stats)
        self.btn_gas_fit.clicked.connect(self.handle_open_gas_fit)
        self.btn_silage_eval.clicked.connect(self.handle_open_silage_eval)
        self.btn_rnaseq.clicked.connect(self.handle_open_rnaseq)

        btn_layout.addWidget(self.btn_16s)
        btn_layout.addWidget(self.btn_feed)
        btn_layout.addWidget(self.btn_stats)
        btn_layout.addWidget(self.btn_gas_fit)
        btn_layout.addWidget(self.btn_silage_eval)
        btn_layout.addWidget(self.btn_rnaseq)

        layout.addStretch()
        layout.addLayout(btn_layout)
        layout.addStretch()

        self.setLayout(layout)

    def handle_open_16s(self):
        if callable(self.on_open_16s):
            self.on_open_16s()

    def handle_open_feed(self):
        if callable(self.on_open_feed):
            self.on_open_feed()

    def handle_open_stats(self):
        if callable(self.on_open_stats):
            self.on_open_stats()

    def handle_open_gas_fit(self):
        if callable(self.on_open_gas_fit):
            self.on_open_gas_fit()

    def handle_open_silage_eval(self):
        if callable(self.on_open_silage_eval):
            self.on_open_silage_eval()

    def handle_open_rnaseq(self):
        if callable(self.on_open_rnaseq):
            self.on_open_rnaseq()

    def apply_style(self):
        self.setStyleSheet("""
            QWidget {
                background-color: #f7f9fc;
                font-family: Microsoft YaHei, Arial;
                font-size: 11pt;
            }

            QLabel#subtitleLabel {
                color: #6b7280;
                font-size: 12pt;
            }

            QPushButton {
                background-color: #4f46e5;
                color: white;
                border: none;
                border-radius: 12px;
                padding: 18px 20px;
                font-size: 14pt;
                font-weight: bold;
                min-height: 56px;
            }

            QPushButton:hover {
                background-color: #4338ca;
            }

            QPushButton:pressed {
                background-color: #3730a3;
            }
        """)
