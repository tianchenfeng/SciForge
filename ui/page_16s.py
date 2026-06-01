import os
import sys
from datetime import datetime

from PySide6.QtWidgets import (
    QWidget, QLabel, QPushButton, QFileDialog,
    QTextEdit, QVBoxLayout, QHBoxLayout, QRadioButton, QCheckBox,
    QMessageBox, QGroupBox, QGridLayout, QFrame, QProgressBar
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from modules.module_16s.config import (
    WINDOW_TITLE,
    MAIN_TITLE,
    SUB_TITLE,
    SCRIPT_DISPLAY_MAP,
    PIPELINE_ORDER,
    DEFAULT_CHECKED,
    FOOTER_LEFT,
    FOOTER_RIGHT,
    NOTICE_TEXT,
    TOOLBOX_VERSION,
    REQUIRED_INPUT_FILES,
    OUTPUT_DIR_HINTS,
    PICRUST2_REQUIRED_FILES,
)
from modules.module_16s.service import SixteenSService
from modules.module_16s.runner import SixteenSRunner


class Page16S(QWidget):
    def __init__(self):
        super().__init__()

        self.input_dir = ""
        self.output_dir = ""

        self.process = None
        self.script_queue = []
        self.current_index = 0
        self.is_running = False

        self.runner = SixteenSRunner(parent=self)

        self.init_ui()
        self.apply_style()

    def get_app_base_dir(self) -> str:
        if getattr(sys, "frozen", False):
            return sys._MEIPASS
        return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(18, 16, 18, 12)
        main_layout.setSpacing(12)

        title_label = QLabel(MAIN_TITLE)
        title_font = QFont("Microsoft YaHei", 16, QFont.Bold)
        title_label.setFont(title_font)
        title_label.setAlignment(Qt.AlignCenter)

        subtitle_label = QLabel(SUB_TITLE)
        subtitle_label.setAlignment(Qt.AlignCenter)
        subtitle_label.setObjectName("subtitleLabel")

        main_layout.addWidget(title_label)
        main_layout.addWidget(subtitle_label)

        path_group = QGroupBox("输入 / 输出设置")
        path_layout = QVBoxLayout()
        path_layout.setSpacing(10)

        input_layout = QHBoxLayout()
        self.input_label = QLabel("输入目录：未选择")
        self.input_label.setObjectName("pathLabel")
        self.btn_input = QPushButton("选择 01_input 文件夹")
        self.btn_input.clicked.connect(self.choose_input_dir)
        input_layout.addWidget(self.input_label, 1)
        input_layout.addWidget(self.btn_input)

        output_layout = QHBoxLayout()
        self.output_label = QLabel("输出目录：未选择")
        self.output_label.setObjectName("pathLabel")
        self.btn_output = QPushButton("选择 02_output 输出位置")
        self.btn_output.clicked.connect(self.choose_output_dir)
        output_layout.addWidget(self.output_label, 1)
        output_layout.addWidget(self.btn_output)

        path_layout.addLayout(input_layout)
        path_layout.addLayout(output_layout)
        path_group.setLayout(path_layout)

        mode_group = QGroupBox("分析模式")
        mode_layout = QHBoxLayout()
        self.radio_single = QRadioButton("单因素分析")
        self.radio_double = QRadioButton("双因素分析")
        self.radio_double.setChecked(True)
        mode_layout.addWidget(self.radio_single)
        mode_layout.addWidget(self.radio_double)
        mode_group.setLayout(mode_layout)

        step_group = QGroupBox("分析步骤")
        step_layout = QVBoxLayout()
        step_layout.setSpacing(8)

        select_btn_layout = QHBoxLayout()
        self.btn_select_all = QPushButton("全选")
        self.btn_select_all.clicked.connect(self.select_all_steps)
        self.btn_select_none = QPushButton("全不选")
        self.btn_select_none.clicked.connect(self.select_no_steps)
        select_btn_layout.addWidget(self.btn_select_all)
        select_btn_layout.addWidget(self.btn_select_none)
        select_btn_layout.addStretch()
        step_layout.addLayout(select_btn_layout)

        self.step_boxes = {}

        step_grid = QGridLayout()
        step_grid.setHorizontalSpacing(28)
        step_grid.setVerticalSpacing(8)

        for i, script_name in enumerate(PIPELINE_ORDER):
            cb = QCheckBox(SCRIPT_DISPLAY_MAP[script_name])
            cb.setChecked(DEFAULT_CHECKED.get(script_name, False))
            self.step_boxes[script_name] = cb
            step_grid.addWidget(cb, i // 2, i % 2)

        step_layout.addLayout(step_grid)
        step_group.setLayout(step_layout)

        btn_layout = QHBoxLayout()
        self.btn_run = QPushButton("开始分析")
        self.btn_run.clicked.connect(self.start_analysis)

        self.btn_cancel = QPushButton("取消分析")
        self.btn_cancel.clicked.connect(self.cancel_analysis)
        self.btn_cancel.setEnabled(False)

        self.btn_open_output = QPushButton("打开输出目录")
        self.btn_open_output.clicked.connect(self.open_output_dir)

        btn_layout.addWidget(self.btn_run)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_open_output)

        progress_group = QGroupBox("运行进度")
        progress_layout = QVBoxLayout()

        self.progress_label = QLabel("当前状态：未开始")
        self.progress_label.setObjectName("progressLabel")

        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat("%p%")

        progress_layout.addWidget(self.progress_label)
        progress_layout.addWidget(self.progress_bar)
        progress_group.setLayout(progress_layout)

        log_group = QGroupBox("运行日志")
        log_layout = QVBoxLayout()

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)

        log_layout.addWidget(self.log_box)
        log_group.setLayout(log_layout)

        # ===== 新增：输出目录结构说明 =====
        output_hint_group = QGroupBox("输出目录说明")
        output_hint_layout = QVBoxLayout()
        output_hint_text = "\n".join([f"{name}    {desc}" for name, desc in OUTPUT_DIR_HINTS])
        self.output_hint_label = QLabel(output_hint_text)
        self.output_hint_label.setObjectName("hintLabel")
        self.output_hint_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        output_hint_layout.addWidget(self.output_hint_label)
        output_hint_group.setLayout(output_hint_layout)

        footer_line = QFrame()
        footer_line.setFrameShape(QFrame.HLine)
        footer_line.setFrameShadow(QFrame.Sunken)

        notice_label = QLabel(NOTICE_TEXT)
        notice_label.setObjectName("noticeLabel")
        notice_label.setAlignment(Qt.AlignCenter)

        footer_layout = QHBoxLayout()
        footer_layout.setContentsMargins(0, 0, 0, 0)

        left_footer = QLabel(FOOTER_LEFT)
        left_footer.setObjectName("footerLabel")

        right_footer = QLabel(FOOTER_RIGHT)
        right_footer.setObjectName("footerLabel")
        right_footer.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

        footer_layout.addWidget(left_footer)
        footer_layout.addStretch()
        footer_layout.addWidget(right_footer)

        main_layout.addWidget(path_group)
        main_layout.addWidget(mode_group)
        main_layout.addWidget(step_group)
        main_layout.addLayout(btn_layout)
        main_layout.addWidget(progress_group)
        main_layout.addWidget(log_group)
        main_layout.addWidget(output_hint_group)
        main_layout.addWidget(footer_line)
        main_layout.addWidget(notice_label)
        main_layout.addLayout(footer_layout)

        self.setLayout(main_layout)

    def apply_style(self):
        self.setStyleSheet("""
            QWidget {
                background-color: #f7f9fc;
                font-family: Microsoft YaHei, Arial;
                font-size: 11pt;
            }

            QGroupBox {
                background-color: white;
                border: 1px solid #d8dee9;
                border-radius: 10px;
                margin-top: 10px;
                padding-top: 12px;
                font-weight: bold;
            }

            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 6px 0 6px;
                color: #2e3440;
            }

            QPushButton {
                background-color: #4f46e5;
                color: white;
                border: none;
                border-radius: 8px;
                padding: 8px 14px;
                font-weight: bold;
                min-height: 22px;
            }

            QPushButton:hover {
                background-color: #4338ca;
            }

            QPushButton:pressed {
                background-color: #3730a3;
            }

            QPushButton:disabled {
                background-color: #a5b4fc;
                color: #eef2ff;
            }

            QTextEdit {
                background-color: #ffffff;
                border: 1px solid #d8dee9;
                border-radius: 8px;
                padding: 6px;
                font-family: Consolas, Courier New, monospace;
                font-size: 10.5pt;
            }

            QLabel#pathLabel {
                color: #374151;
                padding: 6px 8px;
                background-color: #f3f4f6;
                border-radius: 6px;
            }

            QLabel#subtitleLabel {
                color: #6b7280;
                margin-bottom: 4px;
            }

            QLabel#footerLabel {
                color: #6b7280;
                font-size: 10pt;
            }

            QLabel#noticeLabel {
                color: #9ca3af;
                font-size: 10pt;
                padding-top: 2px;
                padding-bottom: 2px;
            }

            QLabel#progressLabel {
                color: #374151;
                font-size: 10.5pt;
            }

            QLabel#hintLabel {
                color: #4b5563;
                background-color: #f9fafb;
                border-radius: 8px;
                padding: 10px;
                font-family: Consolas, Courier New, monospace;
                font-size: 10pt;
            }

            QProgressBar {
                border: 1px solid #d1d5db;
                border-radius: 8px;
                background: #f3f4f6;
                text-align: center;
                height: 22px;
                font-weight: bold;
                color: #111827;
            }

            QProgressBar::chunk {
                background-color: #4f46e5;
                border-radius: 7px;
            }

            QCheckBox {
                spacing: 10px;
                padding: 4px 2px;
                color: #111827;
            }

            QCheckBox::indicator {
                width: 18px;
                height: 18px;
                border: 2px solid #6366f1;
                border-radius: 4px;
                background-color: white;
            }

            QCheckBox::indicator:checked {
                background-color: #4f46e5;
                border: 2px solid #4338ca;
            }

            QCheckBox::indicator:unchecked:hover {
                background-color: #eef2ff;
            }

            QRadioButton {
                spacing: 10px;
                padding: 4px 2px;
                color: #111827;
            }

            QRadioButton::indicator {
                width: 18px;
                height: 18px;
            }

            QRadioButton::indicator:unchecked {
                border: 2px solid #6366f1;
                border-radius: 9px;
                background-color: white;
            }

            QRadioButton::indicator:checked {
                border: 2px solid #4338ca;
                border-radius: 9px;
                background-color: #4f46e5;
            }
        """)

    def choose_input_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "选择 01_input 文件夹")
        if folder:
            self.input_dir = folder
            self.input_label.setText(f"输入目录：{folder}")

    def choose_output_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "选择输出目录")
        if folder:
            self.output_dir = folder
            self.output_label.setText(f"输出目录：{folder}")

    def select_all_steps(self):
        for cb in self.step_boxes.values():
            cb.setChecked(True)

    def select_no_steps(self):
        for cb in self.step_boxes.values():
            cb.setChecked(False)

    def log(self, text: str):
        self.log_box.append(text)

    def set_controls_enabled(self, enabled: bool):
        self.btn_input.setEnabled(enabled)
        self.btn_output.setEnabled(enabled)
        self.radio_single.setEnabled(enabled)
        self.radio_double.setEnabled(enabled)
        self.btn_select_all.setEnabled(enabled)
        self.btn_select_none.setEnabled(enabled)

        for cb in self.step_boxes.values():
            cb.setEnabled(enabled)

        self.btn_run.setEnabled(enabled)
        self.btn_cancel.setEnabled(not enabled)

    def update_progress(self):
        total = len(self.script_queue)
        if total == 0:
            self.progress_bar.setValue(0)
            self.progress_label.setText("当前状态：未开始")
            return

        finished = self.current_index
        percent = int(finished / total * 100)
        self.progress_bar.setValue(percent)

        if finished < total:
            current_script = self.script_queue[finished]
            display_name = SCRIPT_DISPLAY_MAP.get(current_script, current_script)
            self.progress_label.setText(
                f"当前进度：第 {finished + 1} / {total} 步，正在运行：{display_name}"
            )
        else:
            self.progress_label.setText(f"当前进度：第 {total} / {total} 步，已完成")
            self.progress_bar.setValue(100)

    # ===== 新增：输入文件检查 =====
    def validate_input_dir(self) -> tuple[bool, str]:
        if not self.input_dir:
            return False, "请先选择输入目录。"

        qiime2_dir = os.path.join(self.input_dir, "qiime2")

        if not os.path.isdir(qiime2_dir):
            return False, "输入目录缺少 qiime2 子文件夹。"

        missing_files = []
        for file_name in REQUIRED_INPUT_FILES:
            file_path = os.path.join(qiime2_dir, file_name)
            if not os.path.exists(file_path):
                missing_files.append(file_name)

        if missing_files:
            missing_str = "\n".join(missing_files)
            return False, f"qiime2 文件夹中缺少以下必要文件：\n{missing_str}"

        return True, ""

    def start_analysis(self):
        if self.is_running:
            return

        ok, msg = self.validate_input_dir()
        if not ok:
            QMessageBox.warning(self, "提示", msg)
            return

        if not self.output_dir:
            QMessageBox.warning(self, "提示", "请先选择输出目录。")
            return

        mode = "single_factor" if self.radio_single.isChecked() else "two_factor"

        selected_scripts = [
            script_name
            for script_name, cb in self.step_boxes.items()
            if cb.isChecked()
        ]

        if not selected_scripts:
            QMessageBox.warning(self, "提示", "请至少勾选一个分析步骤。")
            return

        # 如果勾选了 05_picrust2.py，则检查 picrust2 文件夹
        if "05_picrust2.py" in selected_scripts:
            picrust2_dir = os.path.join(self.input_dir, "picrust2")

            if not os.path.isdir(picrust2_dir):
                QMessageBox.warning(self, "提示", "输入目录缺少 picrust2 子文件夹。")
                return

            missing_files = []
            for file_name in PICRUST2_REQUIRED_FILES:
                file_path = os.path.join(picrust2_dir, file_name)
                if not os.path.exists(file_path):
                    missing_files.append(file_name)

            if missing_files:
                missing_str = "\n".join(missing_files)
                QMessageBox.warning(self, "提示", f"picrust2 文件夹缺少以下必要文件：\n{missing_str}")
                return

        normalized_scripts = SixteenSService.normalize_selected_scripts(selected_scripts, mode)

        if mode == "single_factor" and "09_diff_visualization.py" in selected_scripts and "09_diff_visualization.py" not in normalized_scripts:
            self.log("提示：单因素模式下，已自动跳过 09 分层差异可视化。")

        self.script_queue = normalized_scripts
        self.current_index = 0
        self.is_running = True

        self.set_controls_enabled(False)
        self.log_box.clear()
        self.progress_bar.setValue(0)

        self.log("===== 16S 分析开始 =====")
        self.log(f"工具箱版本：{TOOLBOX_VERSION}")
        self.log(f"启动时间：{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.log(f"分析模式：{mode}")
        self.log(f"输入目录：{self.input_dir}")
        self.log(f"输出目录：{self.output_dir}")

        self.update_progress()
        self.run_next_script(mode)

    def run_next_script(self, mode: str):
        if self.current_index >= len(self.script_queue):
            self.is_running = False
            self.process = None
            self.update_progress()
            self.log("[DONE] 全部分析完成！")
            self.set_controls_enabled(True)

            log_path = self.save_log_to_file(status="success")
            if log_path:
                self.log(f"[LOG] 日志已保存：{log_path}")

            QMessageBox.information(self, "完成", "全部分析已成功完成！")
            return

        script_name = self.script_queue[self.current_index]
        display_name = SCRIPT_DISPLAY_MAP.get(script_name, script_name)

        self.log(f"\n===== 正在运行：{display_name} =====")
        self.update_progress()

        try:
            self.process = self.runner.create_process(
                script_name=script_name,
                input_dir=self.input_dir,
                output_dir=self.output_dir,
                mode=mode,
                app_base_dir=self.get_app_base_dir(),
            )
        except Exception as e:
            self.is_running = False
            self.process = None
            self.set_controls_enabled(True)
            self.log(f"[ERROR] 无法启动脚本：{display_name}")
            self.log(str(e))
            QMessageBox.critical(self, "错误", f"无法启动脚本：{display_name}\n{e}")
            return

        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(
            lambda exit_code, exit_status: self.handle_finished(exit_code, exit_status, mode)
        )

        self.process.start()

    def handle_stdout(self):
        if self.process is None:
            return
        data = self.process.readAllStandardOutput().data().decode("utf-8", errors="ignore")
        if data.strip():
            self.log(data)

    def handle_stderr(self):
        if self.process is None:
            return
        data = self.process.readAllStandardError().data().decode("utf-8", errors="ignore")
        if data.strip():
            self.log(data)

    def handle_finished(self, exit_code, exit_status, mode: str):
        if self.current_index >= len(self.script_queue):
            return

        script_name = self.script_queue[self.current_index]
        display_name = SCRIPT_DISPLAY_MAP.get(script_name, script_name)

        if exit_code == 0:
            self.log(f"[OK] 完成：{display_name}")
            self.current_index += 1
            self.process = None
            self.update_progress()
            self.run_next_script(mode)
        else:
            self.is_running = False
            self.script_queue = []
            self.process = None
            self.log(f"[ERROR] 失败：{display_name}")
            self.set_controls_enabled(True)
            self.progress_label.setText(f"当前状态：运行失败（{display_name}）")

            log_path = self.save_log_to_file(status="error")
            if log_path:
                self.log(f"[LOG] 错误日志已保存：{log_path}")

            QMessageBox.critical(
                self,
                "错误",
                f"{display_name} 运行失败。\n日志已自动保存，请查看日志窗口或日志文件。"
            )

    def cancel_analysis(self):
        if not self.is_running:
            return

        reply = QMessageBox.question(
            self,
            "确认取消",
            "确定要取消当前分析吗？",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply != QMessageBox.Yes:
            return

        self.log("[WARN] 用户已取消分析。")

        self.is_running = False
        self.script_queue = []

        if self.process is not None:
            self.process.kill()
            self.process.waitForFinished(2000)
            self.process = None

        self.set_controls_enabled(True)
        self.progress_label.setText("当前状态：分析已取消")
        self.progress_bar.setValue(0)

        log_path = self.save_log_to_file(status="cancelled")
        if log_path:
            self.log(f"[LOG] 日志已保存：{log_path}")

    def save_log_to_file(self, status="normal"):
        if not self.output_dir:
            return None

        log_dir = os.path.join(self.output_dir, "_run_logs")
        os.makedirs(log_dir, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = os.path.join(log_dir, f"run_log_{status}_{timestamp}.txt")

        try:
            with open(log_path, "w", encoding="utf-8") as f:
                f.write(self.log_box.toPlainText())
            return log_path
        except Exception:
            return None

    def open_output_dir(self):
        if self.output_dir and os.path.exists(self.output_dir):
            os.startfile(self.output_dir)
        else:
            QMessageBox.warning(self, "提示", "输出目录不存在。")