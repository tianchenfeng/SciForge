from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QFileDialog,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from modules.module_rnaseq_analysis.service import RNASeqAnalysisService


class PageRNASeqAnalysis(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()
        self.apply_style()

    def init_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("RNASeq 综合分析")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Microsoft YaHei", 16, QFont.Bold))
        layout.addWidget(title)

        form_group = QGroupBox("分析设置")
        form = QVBoxLayout()
        self.deg_edit = QLineEdit()
        self.annotation_edit = QLineEdit()
        self.go_edit = QLineEdit()
        self.kegg_edit = QLineEdit()
        self.expression_edit = QLineEdit()
        self.output_edit = QLineEdit()
        self.keyword_edit = QLineEdit()
        self.keyword_edit.setPlaceholderText("例如: cellulose, starch, fermentation")
        self.log2fc_spin = QDoubleSpinBox()
        self.log2fc_spin.setRange(0, 100)
        self.log2fc_spin.setDecimals(3)
        self.log2fc_spin.setValue(1.0)
        self.padj_spin = QDoubleSpinBox()
        self.padj_spin.setRange(0, 1)
        self.padj_spin.setDecimals(4)
        self.padj_spin.setSingleStep(0.01)
        self.padj_spin.setValue(0.05)

        form.addLayout(self._file_row("DEG 表", self.deg_edit, self.choose_deg))
        form.addLayout(self._file_row("Annotation", self.annotation_edit, self.choose_annotation))
        form.addLayout(self._file_row("GO 富集", self.go_edit, self.choose_go))
        form.addLayout(self._file_row("KEGG 富集", self.kegg_edit, self.choose_kegg))
        form.addLayout(self._file_row("表达矩阵", self.expression_edit, self.choose_expression))
        form.addLayout(self._file_row("输出目录", self.output_edit, self.choose_output))

        threshold_row = QHBoxLayout()
        threshold_row.addWidget(QLabel("log2FC 阈值"))
        threshold_row.addWidget(self.log2fc_spin)
        threshold_row.addWidget(QLabel("padj 阈值"))
        threshold_row.addWidget(self.padj_spin)
        threshold_row.addStretch()
        form.addLayout(threshold_row)

        keyword_row = QHBoxLayout()
        keyword_row.addWidget(QLabel("候选关键词"))
        keyword_row.addWidget(self.keyword_edit, 1)
        form.addLayout(keyword_row)

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

    def choose_deg(self):
        self._choose_file(self.deg_edit)

    def choose_annotation(self):
        self._choose_file(self.annotation_edit)

    def choose_go(self):
        self._choose_file(self.go_edit)

    def choose_kegg(self):
        self._choose_file(self.kegg_edit)

    def choose_expression(self):
        self._choose_file(self.expression_edit)

    def _choose_file(self, edit: QLineEdit):
        path, _ = QFileDialog.getOpenFileName(self, "选择输入文件", "", "Data files (*.xlsx *.xls *.csv *.tsv *.txt)")
        if path:
            edit.setText(path)

    def choose_output(self):
        path = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if path:
            self.output_edit.setText(path)

    def log(self, message: str):
        self.log_box.append(message)

    def run_analysis(self):
        try:
            self.log_box.clear()
            result = RNASeqAnalysisService.run(
                self.deg_edit.text().strip(),
                self.output_edit.text().strip(),
                self.log2fc_spin.value(),
                self.padj_spin.value(),
                self.keyword_edit.text().strip(),
                self.annotation_edit.text().strip() or None,
                self.go_edit.text().strip() or None,
                self.kegg_edit.text().strip() or None,
                self.expression_edit.text().strip() or None,
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
            QPushButton { background-color: #7c3aed; color: white; border: none; border-radius: 7px; padding: 8px 14px; font-weight: bold; }
            QPushButton:hover { background-color: #6d28d9; }
            QLineEdit, QTextEdit, QDoubleSpinBox { background-color: white; border: 1px solid #d8dee9; border-radius: 6px; min-height: 30px; }
        """)
