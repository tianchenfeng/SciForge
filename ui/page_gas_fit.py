from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFileDialog,
    QLineEdit, QTextEdit, QComboBox, QGroupBox, QMessageBox
)
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt

from modules.module_gas_fit.service import GasFitService


class PageGasFit(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()
        self.apply_style()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("体外产气拟合")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Microsoft YaHei", 16, QFont.Bold))
        layout.addWidget(title)

        form_group = QGroupBox("分析设置")
        form = QVBoxLayout()

        self.input_edit = QLineEdit()
        self.output_edit = QLineEdit()
        self.model_combo = QComboBox()
        self.model_combo.addItem("指数模型 GP = A*(1-exp(-c*t))", "exponential")
        self.model_combo.addItem("带 Lag 模型", "lag")

        form.addLayout(self._file_row("输入文件", self.input_edit, self.choose_input))
        form.addLayout(self._file_row("输出目录", self.output_edit, self.choose_output))
        model_row = QHBoxLayout()
        model_row.addWidget(QLabel("模型"))
        model_row.addWidget(self.model_combo)
        form.addLayout(model_row)

        self.run_button = QPushButton("开始分析")
        self.run_button.clicked.connect(self.run_analysis)
        form.addWidget(self.run_button)
        form_group.setLayout(form)

        log_group = QGroupBox("日志")
        log_layout = QVBoxLayout()
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        log_layout.addWidget(self.log_box)
        log_group.setLayout(log_layout)

        layout.addWidget(form_group)
        layout.addWidget(log_group, 1)
        self.setLayout(layout)

    def _file_row(self, label: str, edit: QLineEdit, handler):
        row = QHBoxLayout()
        row.addWidget(QLabel(label))
        row.addWidget(edit, 1)
        btn = QPushButton("选择")
        btn.clicked.connect(handler)
        row.addWidget(btn)
        return row

    def choose_input(self):
        path, _ = QFileDialog.getOpenFileName(self, "选择输入文件", "", "Data files (*.xlsx *.xls *.csv)")
        if path:
            self.input_edit.setText(path)

    def choose_output(self):
        path = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if path:
            self.output_edit.setText(path)

    def log(self, message: str):
        self.log_box.append(message)

    def run_analysis(self):
        try:
            self.log_box.clear()
            result = GasFitService.run(
                self.input_edit.text().strip(),
                self.output_edit.text().strip(),
                self.model_combo.currentData(),
                logger=self.log,
            )
            QMessageBox.information(self, "完成", f"分析完成:\n{result['excel']}")
        except Exception as exc:
            QMessageBox.critical(self, "分析失败", str(exc))

    def apply_style(self):
        self.setStyleSheet("""
            QWidget { background-color: #f7f9fc; font-family: Microsoft YaHei, Arial; font-size: 11pt; }
            QGroupBox { background-color: white; border: 1px solid #d8dee9; border-radius: 8px; margin-top: 10px; padding-top: 12px; font-weight: bold; }
            QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 6px; color: #2e3440; }
            QPushButton { background-color: #2563eb; color: white; border: none; border-radius: 7px; padding: 8px 14px; font-weight: bold; }
            QPushButton:hover { background-color: #1d4ed8; }
            QLineEdit, QTextEdit, QComboBox { background-color: white; border: 1px solid #d8dee9; border-radius: 6px; min-height: 30px; }
        """)
