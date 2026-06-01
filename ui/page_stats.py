from pathlib import Path
import os
from datetime import datetime

from PySide6.QtWidgets import (
    QWidget,
    QLabel,
    QPushButton,
    QFileDialog,
    QTextEdit,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QGroupBox,
    QComboBox,
    QLineEdit,
    QMessageBox,
)

MODE_ONE_WAY = "one_way"
MODE_TWO_WAY = "two_way"


class PageStats(QWidget):
    def __init__(self):
        super().__init__()
        self.input_file = ""
        self.output_dir = ""
        self.columns = []
        self.init_ui()
        self.update_mode_ui()

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(18, 16, 18, 16)
        main_layout.setSpacing(12)

        title = QLabel("通用统计分析")
        title.setObjectName("pageTitle")
        main_layout.addWidget(title)

        path_group = QGroupBox("输入与输出")
        path_layout = QGridLayout()
        path_layout.setHorizontalSpacing(10)
        path_layout.setVerticalSpacing(10)

        self.input_edit = QLineEdit()
        self.input_edit.setReadOnly(True)
        self.input_edit.setPlaceholderText("请选择输入表格文件")

        self.output_edit = QLineEdit()
        self.output_edit.setReadOnly(True)
        self.output_edit.setPlaceholderText("请选择结果输出文件夹")

        self.btn_choose_input = QPushButton("选择输入文件")
        self.btn_choose_input.clicked.connect(self.choose_input_file)

        self.btn_choose_output = QPushButton("选择输出目录")
        self.btn_choose_output.clicked.connect(self.choose_output_dir)

        path_layout.addWidget(QLabel("输入文件"), 0, 0)
        path_layout.addWidget(self.input_edit, 0, 1)
        path_layout.addWidget(self.btn_choose_input, 0, 2)
        path_layout.addWidget(QLabel("输出目录"), 1, 0)
        path_layout.addWidget(self.output_edit, 1, 1)
        path_layout.addWidget(self.btn_choose_output, 1, 2)
        path_group.setLayout(path_layout)

        setting_group = QGroupBox("分析设置")
        setting_layout = QGridLayout()
        setting_layout.setHorizontalSpacing(10)
        setting_layout.setVerticalSpacing(10)

        self.combo_mode = QComboBox()
        self.combo_mode.addItem("单因素分析", MODE_ONE_WAY)
        self.combo_mode.addItem("双因素分析", MODE_TWO_WAY)
        self.combo_mode.currentIndexChanged.connect(self.update_mode_ui)

        self.label_factor_a = QLabel("因素列")
        self.label_factor_b = QLabel("因素B")

        self.combo_factor_a = QComboBox()
        self.combo_factor_b = QComboBox()

        self.btn_load_columns = QPushButton("读取表头")
        self.btn_load_columns.clicked.connect(self.load_columns)

        self.metric_hint = QLabel("提示：指标列将根据输入格式自动识别，不需要手动选择。")
        self.metric_hint.setWordWrap(True)

        setting_layout.addWidget(QLabel("分析模式"), 0, 0)
        setting_layout.addWidget(self.combo_mode, 0, 1)
        setting_layout.addWidget(self.btn_load_columns, 0, 2)
        setting_layout.addWidget(self.label_factor_a, 1, 0)
        setting_layout.addWidget(self.combo_factor_a, 1, 1)
        setting_layout.addWidget(self.label_factor_b, 2, 0)
        setting_layout.addWidget(self.combo_factor_b, 2, 1)
        setting_layout.addWidget(self.metric_hint, 3, 0, 1, 3)
        setting_group.setLayout(setting_layout)

        action_layout = QHBoxLayout()
        action_layout.addStretch()
        self.btn_run = QPushButton("开始分析")
        self.btn_run.clicked.connect(self.run_analysis)
        action_layout.addWidget(self.btn_run)

        log_group = QGroupBox()
        log_layout = QVBoxLayout()

        log_title_layout = QHBoxLayout()
        log_title = QLabel("日志输出")
        self.btn_clear_log = QPushButton("清空日志")
        self.btn_clear_log.clicked.connect(self.log_box_clear)

        log_title_layout.addWidget(log_title)
        log_title_layout.addStretch()
        log_title_layout.addWidget(self.btn_clear_log)

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setFixedHeight(260)

        log_layout.addLayout(log_title_layout)
        log_layout.addWidget(self.log_box)
        log_group.setLayout(log_layout)

        main_layout.addWidget(path_group)
        main_layout.addWidget(setting_group)
        main_layout.addLayout(action_layout)
        main_layout.addWidget(log_group)

        self.setLayout(main_layout)

    def current_mode(self):
        return self.combo_mode.currentData()

    def update_mode_ui(self):
        is_one_way = self.current_mode() == MODE_ONE_WAY

        if is_one_way:
            self.label_factor_a.setText("分组列")
            self.label_factor_b.setVisible(False)
            self.combo_factor_b.setVisible(False)
            self.metric_hint.setText("单因素格式：第1列为分组列，第2列及以后为指标列。")
        else:
            self.label_factor_a.setText("因素A")
            self.label_factor_b.setText("因素B")
            self.label_factor_b.setVisible(True)
            self.combo_factor_b.setVisible(True)
            self.metric_hint.setText("双因素格式：第1列为因素A，第2列为因素B，第3列及以后为指标列。")

    def choose_input_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "选择输入文件",
            "",
            "表格文件 (*.xlsx *.csv *.tsv *.txt);;所有文件 (*.*)",
        )
        if not file_path:
            return
        self.input_file = file_path
        self.input_edit.setText(file_path)
        self.log(f"已选择输入文件：{file_path}")

    def choose_output_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "选择输出目录", "")
        if not folder:
            return
        self.output_dir = folder
        self.output_edit.setText(folder)
        self.log(f"已选择输出目录：{folder}")

    def load_columns(self):
        if not self.input_file:
            self.warn("请先选择输入文件。")
            return

        try:
            from modules.module_stats.service import get_table_columns
            self.columns = get_table_columns(self.input_file)
        except Exception as exc:
            self.log(f"读取表头失败：{exc}")
            self.warn(f"读取表头失败：{exc}")
            return

        self.set_columns(self.columns)
        self.log(f"表头读取成功，共检测到 {len(self.columns)} 列。")
        self.log_detected_structure()

    def set_columns(self, columns):
        clean_columns = [str(col) for col in columns]

        for combo in (self.combo_factor_a, self.combo_factor_b):
            combo.clear()
            combo.addItems(clean_columns)

        if self.combo_factor_a.count() >= 1:
            self.combo_factor_a.setCurrentIndex(0)

        if self.combo_factor_b.count() >= 2:
            self.combo_factor_b.setCurrentIndex(1)
        elif self.combo_factor_b.count() >= 1:
            self.combo_factor_b.setCurrentIndex(0)

    def log_detected_structure(self):
        if not self.columns:
            return

        if self.current_mode() == MODE_ONE_WAY:
            factor = self.columns[0] if len(self.columns) >= 1 else ""
            metrics = self.columns[1:] if len(self.columns) >= 2 else []

            self.log("当前识别规则：单因素分析")
            self.log(f"分组列：{factor}")
            self.log(f"候选指标列数量：{len(metrics)}")
            self.log("候选指标列：" + "，".join(map(str, metrics)) if metrics else "未检测到候选指标列，请检查表格格式。")
        else:
            factor_a = self.columns[0] if len(self.columns) >= 1 else ""
            factor_b = self.columns[1] if len(self.columns) >= 2 else ""
            metrics = self.columns[2:] if len(self.columns) >= 3 else []

            self.log("当前识别规则：双因素分析")
            self.log(f"因素A：{factor_a}")
            self.log(f"因素B：{factor_b}")
            self.log(f"候选指标列数量：{len(metrics)}")
            self.log("候选指标列：" + "，".join(map(str, metrics)) if metrics else "未检测到候选指标列，请检查表格格式。")

    def run_analysis(self):
        if not self.input_file:
            self.warn("请先选择输入文件。")
            return

        if not self.output_dir:
            self.warn("请先选择输出目录。")
            return

        if not self.columns:
            self.warn("请先点击“读取表头”。")
            return

        self.btn_run.setEnabled(False)
        self.log("开始批量统计分析...")

        try:
            if self.current_mode() == MODE_ONE_WAY:
                from modules.module_stats.service import run_one_way_batch_analysis
                result = run_one_way_batch_analysis(
                    input_file=self.input_file,
                    output_dir=self.output_dir,
                    factor=None,
                )
            else:
                from modules.module_stats.service import run_two_way_batch_analysis
                result = run_two_way_batch_analysis(
                    input_file=self.input_file,
                    output_dir=self.output_dir,
                    factor1=None,
                    factor2=None,
                )
        except Exception as exc:
            self.btn_run.setEnabled(True)
            self.log(f"统计分析失败：{exc}")
            self.warn(f"统计分析失败：{exc}")
            return

        self.btn_run.setEnabled(True)
        self.handle_result(result)

    def handle_result(self, result):
        if not isinstance(result, dict):
            self.log(f"分析结果异常：{result}")
            self.warn("分析结果异常。")
            return

        message = result.get("message", "")
        if message:
            if result.get("analysis_type") == "one_way":
                self.log("单因素批量统计分析完成。")
            elif result.get("analysis_type") == "two_way":
                self.log("双因素批量统计分析完成。")
            else:
                self.log(message)

        if result.get("success"):
            metrics = result.get("metrics", [])
            word_path = result.get("word_path", "")
            excel_path = result.get("excel_path", "")

            self.log(f"检测到指标数量：{len(metrics)}")

            if metrics:
                self.log("指标名称：" + "，".join(map(str, metrics)))

            if word_path:
                self.log(f"Word报告：{Path(word_path)}")

            if excel_path:
                self.log(f"Excel结果：{Path(excel_path)}")

            QMessageBox.information(
                self,
                "完成",
                "统计分析完成。\n\nWord报告和Excel结果已生成。"
            )

            self.open_output_dir()
        else:
            self.warn(message or "统计分析失败。")

    def open_output_dir(self):
        if not self.output_dir:
            return

        try:
            os.startfile(self.output_dir)
        except Exception:
            pass

    def log(self, message):
        now = datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{now}] {message}")
        scrollbar = self.log_box.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    def log_box_clear(self):
        self.log_box.clear()

    def warn(self, message):
        QMessageBox.warning(self, "提示", str(message))