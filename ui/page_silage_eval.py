from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFileDialog,
    QLineEdit, QTextEdit, QGroupBox, QMessageBox
)
from PySide6.QtGui import QFont
from PySide6.QtCore import Qt

from modules.module_silage_eval.service import SilageEvalService


class PageSilageEval(QWidget):
    def __init__(self):
        super().__init__()
        self.directions = {}
        self.init_ui()
        self.apply_style()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("青贮综合评价")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Microsoft YaHei", 16, QFont.Bold))
        layout.addWidget(title)

        form_group = QGroupBox("分析设置")
        form = QVBoxLayout()
        self.input_edit = QLineEdit()
        self.output_edit = QLineEdit()
        form.addLayout(self._file_row("输入文件", self.input_edit, self.choose_input))
        form.addLayout(self._file_row("输出目录", self.output_edit, self.choose_output))
        self.inspect_button = QPushButton("识别指标方向")
        self.inspect_button.clicked.connect(self.inspect_directions)
        self.run_button = QPushButton("开始分析")
        self.run_button.clicked.connect(self.run_analysis)
        btn_row = QHBoxLayout()
        btn_row.addWidget(self.inspect_button)
        btn_row.addWidget(self.run_button)
        btn_row.addStretch()
        form.addLayout(btn_row)
        form_group.setLayout(form)

        self.direction_box = QTextEdit()
        self.direction_box.setReadOnly(True)
        direction_group = QGroupBox("识别到的指标方向")
        direction_layout = QVBoxLayout()
        direction_layout.addWidget(self.direction_box)
        direction_group.setLayout(direction_layout)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        log_group = QGroupBox("日志")
        log_layout = QVBoxLayout()
        log_layout.addWidget(self.log_box)
        log_group.setLayout(log_layout)

        layout.addWidget(form_group)
        layout.addWidget(direction_group)
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
            self.inspect_directions()

    def choose_output(self):
        path = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if path:
            self.output_edit.setText(path)

    def inspect_directions(self):
        try:
            self.directions = SilageEvalService.inspect_directions(self.input_edit.text().strip())
            lines = [f"{key}: {'负向' if value == 'negative' else '正向'}" for key, value in self.directions.items()]
            self.direction_box.setPlainText("\n".join(lines))
        except Exception as exc:
            QMessageBox.warning(self, "识别失败", str(exc))

    def log(self, message: str):
        self.log_box.append(message)

    def run_analysis(self):
        try:
            self.log_box.clear()
            result = SilageEvalService.run(
                self.input_edit.text().strip(),
                self.output_edit.text().strip(),
                self.directions or None,
                logger=self.log,
            )
            QMessageBox.information(self, "完成", f"分析完成:\n{result['excel']}\n{result['word']}")
        except Exception as exc:
            QMessageBox.critical(self, "分析失败", str(exc))

    def apply_style(self):
        self.setStyleSheet("""
            QWidget { background-color: #f7f9fc; font-family: Microsoft YaHei, Arial; font-size: 11pt; }
            QGroupBox { background-color: white; border: 1px solid #d8dee9; border-radius: 8px; margin-top: 10px; padding-top: 12px; font-weight: bold; }
            QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 6px; color: #2e3440; }
            QPushButton { background-color: #0f766e; color: white; border: none; border-radius: 7px; padding: 8px 14px; font-weight: bold; }
            QPushButton:hover { background-color: #0d665f; }
            QLineEdit, QTextEdit { background-color: white; border: 1px solid #d8dee9; border-radius: 6px; min-height: 30px; }
        """)
