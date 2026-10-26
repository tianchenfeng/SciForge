from __future__ import annotations

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QTableWidget,
    QTableWidgetItem, QHeaderView, QGroupBox, QGridLayout, QCheckBox,
    QDoubleSpinBox, QMessageBox, QTextEdit, QComboBox, QAbstractSpinBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from modules.module_feed_formula.service import FeedFormulaService


class PageFeedFormula(QWidget):
    COLUMN_KEYS = ["name", "type", "cost", "dm", "cp", "starch", "ndf", "adf", "ca", "p"]
    COLUMN_LABELS = ["原料名称", "类型", "价格", "DM", "CP", "淀粉", "NDF", "ADF", "Ca", "P"]

    DEFAULT_ROWS = [
        {"name": "原料A", "type": "roughage", "cost": 2.5, "dm": 90, "cp": 18, "starch": 1.2, "ndf": 42, "adf": 30, "ca": 1.4, "p": 0.25},
        {"name": "原料B", "type": "concentrate", "cost": 3.2, "dm": 89, "cp": 46.5, "starch": 5.0, "ndf": 12, "adf": 5, "ca": 0.35, "p": 0.65},
        {"name": "原料C", "type": "concentrate", "cost": 2.0, "dm": 88, "cp": 8.5, "starch": 70, "ndf": 10, "adf": 3, "ca": 0.03, "p": 0.28},
    ]

    def __init__(self):
        super().__init__()
        self.constraint_widgets = {}
        self.init_ui()
        self.apply_style()
        self.load_default_rows()

    def init_ui(self):
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(12)

        title = QLabel("智能饲料配方")
        title.setFont(QFont("Microsoft YaHei", 16, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)

        subtitle = QLabel("自由配方版：最低成本 + 营养约束 + 粗料比例约束")
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setObjectName("subtitleLabel")

        main_layout.addWidget(title)
        main_layout.addWidget(subtitle)

        # 原料表
        table_group = QGroupBox("1. 原料成分录入（干物质基础 %）")
        table_layout = QVBoxLayout()

        btn_row = QHBoxLayout()
        self.btn_add_row = QPushButton("添加原料")
        self.btn_add_row.clicked.connect(self.add_empty_row)

        self.btn_delete_row = QPushButton("删除选中行")
        self.btn_delete_row.clicked.connect(self.delete_selected_rows)

        self.btn_reset_demo = QPushButton("载入示例")
        self.btn_reset_demo.clicked.connect(self.reset_demo)

        btn_row.addWidget(self.btn_add_row)
        btn_row.addWidget(self.btn_delete_row)
        btn_row.addWidget(self.btn_reset_demo)
        btn_row.addStretch()

        self.table = QTableWidget(0, len(self.COLUMN_LABELS))
        self.table.setHorizontalHeaderLabels(self.COLUMN_LABELS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setMinimumSectionSize(82)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Fixed)
        self.table.setColumnWidth(1, 150)
        self.table.verticalHeader().setVisible(True)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)

        table_layout.addLayout(btn_row)
        table_layout.addWidget(self.table)
        table_group.setLayout(table_layout)

        # 约束区
        constraint_group = QGroupBox("2. 约束条件设置（Target ± Error）")
        constraint_layout = QGridLayout()
        self.constraint_layout = constraint_layout
        self.constraint_boxes = []
        self.constraint_columns = 2
        constraint_layout.setHorizontalSpacing(18)
        constraint_layout.setVerticalSpacing(10)

        constraint_defs = [
            ("cp", "粗蛋白 CP", True, 17.5, 1.0, 0.1, 0.1),
            ("starch", "淀粉", True, 24.0, 2.0, 0.1, 0.1),
            ("ndf", "NDF", False, 30.0, 2.0, 0.1, 0.1),
            ("adf", "ADF", False, 19.0, 2.0, 0.1, 0.1),
            ("ca", "Ca", False, 0.6, 0.1, 0.01, 0.01),
            ("p", "P", False, 0.4, 0.05, 0.01, 0.01),
            ("ratio", "粗料占比", True, 50.0, 5.0, 1.0, 0.5),
        ]

        for idx, (key, label_text, default_enable, default_target, default_error, target_step, error_step) in enumerate(constraint_defs):
            box = self._build_constraint_box(
                key, label_text, default_enable, default_target, default_error, target_step, error_step
            )
            constraint_layout.addWidget(box, idx // 2, idx % 2)
            self.constraint_boxes.append(box)

        constraint_group.setLayout(constraint_layout)

        # 运行按钮
        action_layout = QHBoxLayout()
        self.btn_calculate = QPushButton("开始计算")
        self.btn_calculate.clicked.connect(self.calculate_formula)

        action_layout.addWidget(self.btn_calculate)
        action_layout.addStretch()

        # 结果区
        result_group = QGroupBox("3. 计算结果")
        result_layout = QVBoxLayout()

        self.result_box = QTextEdit()
        self.result_box.setReadOnly(True)

        result_layout.addWidget(self.result_box)
        result_group.setLayout(result_layout)

        main_layout.addWidget(table_group)
        main_layout.addWidget(constraint_group)
        main_layout.addLayout(action_layout)
        main_layout.addWidget(result_group)

        self.setLayout(main_layout)

    def _build_constraint_box(
            self,
            key: str,
            label_text: str,
            default_enable: bool,
            default_target: float,
            default_error: float,
            target_step: float,
            error_step: float,
    ) -> QGroupBox:
        box = QGroupBox(label_text)
        layout = QGridLayout()
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setHorizontalSpacing(12)
        layout.setVerticalSpacing(8)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 1)

        enable_cb = QCheckBox("启用")
        enable_cb.setChecked(default_enable)

        target_label = QLabel("目标值")
        target_label.setMinimumWidth(52)
        target_label.setAttribute(Qt.WA_TransparentForMouseEvents, True)

        target_spin = QDoubleSpinBox()
        target_spin.setRange(0, 999999)
        target_spin.setDecimals(4)
        target_spin.setSingleStep(target_step)
        target_spin.setValue(default_target)
        target_spin.setKeyboardTracking(False)
        target_spin.setMinimumHeight(32)
        target_spin.setMinimumWidth(105)
        target_spin.setAlignment(Qt.AlignLeft)
        target_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        from PySide6.QtGui import QFont

        def _style_pm_button(btn: QPushButton):
            btn.setFixedSize(32, 32)
            btn.setProperty("compact", True)
            btn.setFont(QFont("Microsoft YaHei", 16, QFont.Bold))  # 放大符号
            btn.setStyleSheet("""
                QPushButton {
                    background-color: #4f46e5;
                    color: white;
                    border: none;
                    border-radius: 10px;
                }
                QPushButton:hover {
                    background-color: #4338ca;
                }
                QPushButton:pressed {
                    background-color: #3730a3;
                }
            """)

        target_minus = QPushButton("−")  # 注意：用真正的减号（更居中）
        target_plus = QPushButton("+")
        _style_pm_button(target_minus)
        _style_pm_button(target_plus)

        error_minus = QPushButton("−")
        error_plus = QPushButton("+")
        _style_pm_button(error_minus)
        _style_pm_button(error_plus)

        target_minus.clicked.connect(target_spin.stepDown)
        target_plus.clicked.connect(target_spin.stepUp)

        target_row = QHBoxLayout()
        target_row.setContentsMargins(0, 0, 0, 0)
        target_row.setSpacing(6)
        target_row.addWidget(target_spin)
        target_row.addWidget(target_minus)
        target_row.addWidget(target_plus)

        error_label = QLabel("误差")
        error_label.setMinimumWidth(52)
        error_label.setAttribute(Qt.WA_TransparentForMouseEvents, True)

        error_spin = QDoubleSpinBox()
        error_spin.setRange(0, 999999)
        error_spin.setDecimals(4)
        error_spin.setSingleStep(error_step)
        error_spin.setValue(default_error)
        error_spin.setKeyboardTracking(False)
        error_spin.setMinimumHeight(32)
        error_spin.setMinimumWidth(105)
        error_spin.setAlignment(Qt.AlignLeft)
        error_spin.setButtonSymbols(QAbstractSpinBox.NoButtons)

        error_minus.clicked.connect(error_spin.stepDown)
        error_plus.clicked.connect(error_spin.stepUp)

        error_row = QHBoxLayout()
        error_row.setContentsMargins(0, 0, 0, 0)
        error_row.setSpacing(6)
        error_row.addWidget(error_spin)
        error_row.addWidget(error_minus)
        error_row.addWidget(error_plus)

        layout.addWidget(enable_cb, 0, 0, 1, 2)
        layout.addWidget(target_label, 1, 0)
        layout.addLayout(target_row, 1, 1)
        layout.addWidget(error_label, 2, 0)
        layout.addLayout(error_row, 2, 1)

        box.setLayout(layout)

        self.constraint_widgets[key] = {
            "enable": enable_cb,
            "target": target_spin,
            "error": error_spin,
            "target_minus": target_minus,
            "target_plus": target_plus,
            "error_minus": error_minus,
            "error_plus": error_plus,
        }

        enable_cb.toggled.connect(
            lambda checked,
                   t=target_spin,
                   e=error_spin,
                   tm=target_minus,
                   tp=target_plus,
                   em=error_minus,
                   ep=error_plus: self._toggle_constraint_inputs(checked, t, e, tm, tp, em, ep)
        )
        self._toggle_constraint_inputs(
            default_enable,
            target_spin,
            error_spin,
            target_minus,
            target_plus,
            error_minus,
            error_plus,
        )

        return box

    def resizeEvent(self, event):
        columns = 2 if event.size().width() >= 850 else 1
        if columns != self.constraint_columns:
            self.constraint_columns = columns
            for box in self.constraint_boxes:
                self.constraint_layout.removeWidget(box)
            for index, box in enumerate(self.constraint_boxes):
                self.constraint_layout.addWidget(box, index // columns, index % columns)
        super().resizeEvent(event)

    def _toggle_constraint_inputs(
            self,
            checked: bool,
            target_spin: QDoubleSpinBox,
            error_spin: QDoubleSpinBox,
            target_minus: QPushButton,
            target_plus: QPushButton,
            error_minus: QPushButton,
            error_plus: QPushButton,
    ):
        target_spin.setEnabled(checked)
        error_spin.setEnabled(checked)
        target_minus.setEnabled(checked)
        target_plus.setEnabled(checked)
        error_minus.setEnabled(checked)
        error_plus.setEnabled(checked)

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

            QTextEdit, QTableWidget, QComboBox {
                background-color: white;
                border: 1px solid #d8dee9;
                border-radius: 6px;
            }

            QDoubleSpinBox {
                background-color: white;
                border: 1px solid #d8dee9;
                border-radius: 6px;
                min-height: 30px;
                min-width: 160px;
                padding-left: 6px;
            }

            QLabel#subtitleLabel {
                color: #6b7280;
            }
        """)

    def load_default_rows(self):
        self.table.setRowCount(0)
        for row_data in self.DEFAULT_ROWS:
            self.add_row(row_data)

    def reset_demo(self):
        self.load_default_rows()
        self.result_box.clear()

    def add_empty_row(self):
        row_data = {
            "name": "新原料",
            "type": "concentrate",
            "cost": 0,
            "dm": 88,
            "cp": 0,
            "starch": 0,
            "ndf": 0,
            "adf": 0,
            "ca": 0,
            "p": 0,
        }
        self.add_row(row_data)

    def add_row(self, row_data: dict):
        row = self.table.rowCount()
        self.table.insertRow(row)

        self.table.setItem(row, 0, QTableWidgetItem(str(row_data["name"])))

        combo = QComboBox()
        combo.addItems(["concentrate", "roughage"])
        combo.setCurrentText(str(row_data["type"]))
        self.table.setCellWidget(row, 1, combo)

        for col, key in enumerate(self.COLUMN_KEYS[2:], start=2):
            item = QTableWidgetItem(str(row_data[key]))
            self.table.setItem(row, col, item)

    def delete_selected_rows(self):
        selected_rows = sorted({idx.row() for idx in self.table.selectedIndexes()}, reverse=True)
        for row in selected_rows:
            self.table.removeRow(row)

    def collect_ingredients(self):
        ingredients = []
        for row in range(self.table.rowCount()):
            name_item = self.table.item(row, 0)
            name = name_item.text().strip() if name_item else ""

            type_widget = self.table.cellWidget(row, 1)
            type_value = type_widget.currentText() if type_widget else "concentrate"

            row_data = {
                "name": name,
                "type": type_value,
            }

            for col, key in enumerate(self.COLUMN_KEYS[2:], start=2):
                item = self.table.item(row, col)
                text = item.text().strip() if item else "0"
                row_data[key] = float(text) if text else 0.0

            ingredients.append(row_data)

        return ingredients

    def collect_constraints(self):
        constraints = {}
        for key, widgets in self.constraint_widgets.items():
            constraints[key] = {
                "enable": widgets["enable"].isChecked(),
                "target": widgets["target"].value(),
                "error": widgets["error"].value(),
            }
        return constraints

    def calculate_formula(self):
        try:
            ingredients = self.collect_ingredients()
            constraints = self.collect_constraints()

            result = FeedFormulaService.run_formula(ingredients, constraints)

            if not result["success"]:
                QMessageBox.warning(self, "计算失败", result["message"])
                return

            self.show_result(result)

        except ValueError as e:
            QMessageBox.warning(self, "输入错误", str(e))
        except Exception as e:
            QMessageBox.critical(self, "程序错误", str(e))

    def show_result(self, result: dict):
        mix = result["mix"]
        stats = result["stats"]

        lines = []
        lines.append("===== 配方结果 =====\n")
        for row in mix:
            lines.append(f"{row['name']:<12} {row['amount']:>8.2f} %")

        lines.append("\n===== 营养验算 =====\n")
        lines.append(f"总成本      : {stats['cost']:.4f}")
        lines.append(f"粗料占比    : {stats['roughage']:.4f} %")
        lines.append(f"粗蛋白 CP   : {stats['cp']:.4f} %")
        lines.append(f"淀粉        : {stats['starch']:.4f} %")
        lines.append(f"NDF         : {stats['ndf']:.4f} %")
        lines.append(f"ADF         : {stats['adf']:.4f} %")
        lines.append(f"Ca          : {stats['ca']:.4f} %")
        lines.append(f"P           : {stats['p']:.4f} %")

        lines.append(f"\n状态：{result['solver_status']}")

        self.result_box.setPlainText("\n".join(lines))
