from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import tempfile

from PySide6.QtCore import QProcess, QUrl, Qt
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar, QPushButton,
    QScrollArea, QSpinBox, QSplitter, QTabWidget, QTableWidget, QTableWidgetItem, QTextBrowser,
    QVBoxLayout, QWidget,
)


class PageMetabolomics(QWidget):
    """SciForge metabolomics workbench with a process-isolated analysis runner."""

    RESULT_AREAS = [
        "数据", "QC", "预处理", "PCA", "PLS-DA/OPLS-DA", "单因素", "双因素",
        "简单效应", "差异代谢物", "热图", "注释", "通路", "相关性", "导出报告",
    ]

    def __init__(self):
        super().__init__()
        self.process: QProcess | None = None
        self.config_path: str | None = None
        self.result_dir: str | None = None
        self.process_buffer = ""
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout(self)
        title = QLabel("代谢组综合分析")
        title.setStyleSheet("font-size: 18pt; font-weight: 700; color: #111827;")
        root.addWidget(title)

        splitter = QSplitter()
        self.splitter = splitter
        splitter.setChildrenCollapsible(False)
        controls = QWidget(); controls_layout = QVBoxLayout(controls)

        io_box = QGroupBox("1. 数据与输出")
        io_form = QFormLayout(io_box)
        self.input_edit = QLineEdit(); self.output_edit = QLineEdit()
        input_row = QHBoxLayout(); input_row.addWidget(self.input_edit)
        choose_dir = QPushButton("选择目录"); choose_file = QPushButton("选择矩阵")
        choose_dir.clicked.connect(self._choose_input_dir); choose_file.clicked.connect(self._choose_input_file)
        input_row.addWidget(choose_dir); input_row.addWidget(choose_file)
        output_row = QHBoxLayout(); output_row.addWidget(self.output_edit)
        choose_output = QPushButton("选择输出目录"); choose_output.clicked.connect(self._choose_output)
        output_row.addWidget(choose_output)
        io_form.addRow("输入", input_row); io_form.addRow("输出", output_row)
        self.phenotype_edit = QLineEdit(); phenotype_row = QHBoxLayout(); phenotype_row.addWidget(self.phenotype_edit)
        phenotype_button = QPushButton("选择"); phenotype_button.clicked.connect(lambda: self._choose_optional_file(self.phenotype_edit, "选择表型表")); phenotype_row.addWidget(phenotype_button)
        self.taxon_edit = QLineEdit(); taxon_row = QHBoxLayout(); taxon_row.addWidget(self.taxon_edit)
        taxon_button = QPushButton("选择"); taxon_button.clicked.connect(lambda: self._choose_optional_file(self.taxon_edit, "选择 16S abundance 表")); taxon_row.addWidget(taxon_button)
        io_form.addRow("表型表（可选）", phenotype_row); io_form.addRow("16S 表（可选）", taxon_row)
        self.inspect_button = QPushButton("识别并校验数据"); self.inspect_button.clicked.connect(self.inspect_data)
        io_form.addRow(self.inspect_button)
        controls_layout.addWidget(io_box)

        design_box = QGroupBox("2. 设计与预处理")
        design = QFormLayout(design_box)
        self.mode_combo = QComboBox(); self.mode_combo.addItem("单因素", "one_factor"); self.mode_combo.addItem("双因素（含交互与简单效应）", "two_factor")
        self.factor_a = QLineEdit("Group"); self.factor_b = QLineEdit("Factor B")
        self.scale_combo = QComboBox(); self.scale_combo.addItem("请选择（必须确认）", "unknown"); self.scale_combo.addItem("原始丰度", "raw"); self.scale_combo.addItem("log2 丰度", "log2"); self.scale_combo.addItem("log10 丰度", "log10")
        self.missing_combo = QComboBox()
        for label, value in [("不填补", "none"), ("全局最小正值的一半", "half_min"), ("组内最小正值的一半", "group_half_min"), ("中位数", "median"), ("KNN", "knn")]: self.missing_combo.addItem(label, value)
        self.transform_combo = QComboBox()
        for value in ["none", "log2", "log10"]: self.transform_combo.addItem(value, value)
        self.scaling_combo = QComboBox()
        for label, value in [("Pareto", "pareto"), ("UV / autoscaling", "uv"), ("中心化", "center"), ("不缩放", "none")]: self.scaling_combo.addItem(label, value)
        design.addRow("分析设计", self.mode_combo); design.addRow("因子 A 列", self.factor_a); design.addRow("因子 B 列", self.factor_b)
        design.addRow("输入尺度", self.scale_combo); design.addRow("缺失值", self.missing_combo); design.addRow("变换", self.transform_combo); design.addRow("缩放", self.scaling_combo)
        controls_layout.addWidget(design_box)

        stats_box = QGroupBox("3. 统计与验证")
        stats = QFormLayout(stats_box)
        self.method_combo = QComboBox(); self.method_combo.addItem("ANOVA", "anova"); self.method_combo.addItem("Welch ANOVA", "welch"); self.method_combo.addItem("Kruskal-Wallis", "kruskal")
        self.posthoc_combo = QComboBox(); self.posthoc_combo.addItem("Tukey HSD", "tukey"); self.posthoc_combo.addItem("Games-Howell", "games_howell"); self.posthoc_combo.addItem("Dunn", "dunn")
        self.fdr_spin = QDoubleSpinBox(); self.fdr_spin.setRange(.001, .5); self.fdr_spin.setDecimals(3); self.fdr_spin.setValue(.05)
        self.fc_spin = QDoubleSpinBox(); self.fc_spin.setRange(1, 100); self.fc_spin.setValue(1.2)
        self.pls_check = QCheckBox("运行 PLS-DA（嵌套交叉验证）"); self.pls_check.setChecked(True)
        self.permutation_spin = QSpinBox(); self.permutation_spin.setRange(0, 10000); self.permutation_spin.setValue(200)
        self.component_spin = QSpinBox(); self.component_spin.setRange(1, 20); self.component_spin.setValue(5)
        self.cv_spin = QSpinBox(); self.cv_spin.setRange(2, 20); self.cv_spin.setValue(5)
        self.heatmap_spin = QSpinBox(); self.heatmap_spin.setRange(5, 200); self.heatmap_spin.setValue(30)
        self.seed_spin = QSpinBox(); self.seed_spin.setRange(0, 2147483647); self.seed_spin.setValue(2026)
        stats.addRow("总体检验", self.method_combo); stats.addRow("事后比较", self.posthoc_combo); stats.addRow("FDR 阈值", self.fdr_spin); stats.addRow("FC 阈值", self.fc_spin)
        stats.addRow(self.pls_check); stats.addRow("最大 PLS 组件", self.component_spin); stats.addRow("CV 折数", self.cv_spin); stats.addRow("置换次数", self.permutation_spin); stats.addRow("热图 top N", self.heatmap_spin); stats.addRow("随机种子", self.seed_spin)
        controls_layout.addWidget(stats_box)

        run_row = QHBoxLayout()
        self.run_button = QPushButton("开始分析"); self.cancel_button = QPushButton("取消"); self.cancel_button.setEnabled(False)
        self.run_button.clicked.connect(self.start_analysis); self.cancel_button.clicked.connect(self.cancel_analysis)
        run_row.addWidget(self.run_button); run_row.addWidget(self.cancel_button); controls_layout.addLayout(run_row)
        self.progress = QProgressBar(); controls_layout.addWidget(self.progress)
        self.stage_label = QLabel("等待数据校验"); controls_layout.addWidget(self.stage_label); controls_layout.addStretch()

        right = QTabWidget()
        metadata_page = QWidget(); metadata_layout = QVBoxLayout(metadata_page)
        self.metadata_table = QTableWidget(0, 5); self.metadata_table.setHorizontalHeaderLabels(["Sample ID", "Sample type", "Group", "Factor A", "Factor B"])
        metadata_layout.addWidget(QLabel("识别后可直接修改分组与因素；运行时会冻结到 manifest。")); metadata_layout.addWidget(self.metadata_table)
        right.addTab(metadata_page, "样本设计")
        self.log = QTextBrowser(); right.addTab(self.log, "运行日志")
        self.results_tabs = QTabWidget()
        self.result_browsers = []
        for area in self.RESULT_AREAS:
            browser = QTextBrowser(); browser.setHtml(f"<h3>{area}</h3><p>完成分析后，请在结果目录及综合报告中查看。</p>")
            self.result_browsers.append(browser)
            self.results_tabs.addTab(browser, area)
        right.addTab(self.results_tabs, "结果浏览")
        self.open_result_button = QPushButton("打开结果目录"); self.open_result_button.setEnabled(False); self.open_result_button.clicked.connect(self.open_result)
        metadata_layout.addWidget(self.open_result_button)

        controls_scroll = QScrollArea()
        controls_scroll.setObjectName("controlsScroll")
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setWidget(controls)
        controls_scroll.setMinimumSize(420, 360)
        right.setMinimumSize(340, 320)
        splitter.addWidget(controls_scroll); splitter.addWidget(right); splitter.setSizes([470, 650])
        root.addWidget(splitter)

    def resizeEvent(self, event):
        orientation = Qt.Horizontal if event.size().width() >= 1050 else Qt.Vertical
        if self.splitter.orientation() != orientation:
            self.splitter.setOrientation(orientation)
        super().resizeEvent(event)

    def _choose_input_dir(self):
        value = QFileDialog.getExistingDirectory(self, "选择 workflow_results 目录")
        if value: self.input_edit.setText(value)

    def _choose_input_file(self):
        value, _ = QFileDialog.getOpenFileName(self, "选择丰度矩阵", filter="Data (*.csv *.tsv *.txt *.xlsx *.xls)")
        if value: self.input_edit.setText(value)

    def _choose_output(self):
        value = QFileDialog.getExistingDirectory(self, "选择输出父目录")
        if value: self.output_edit.setText(value)

    def _choose_optional_file(self, target, title):
        value, _ = QFileDialog.getOpenFileName(self, title, filter="Data (*.csv *.tsv *.txt *.xlsx *.xls)")
        if value: target.setText(value)

    def inspect_data(self):
        source = self.input_edit.text().strip()
        if not source or not Path(source).exists():
            QMessageBox.warning(self, "输入无效", "请选择存在的 workflow_results 目录或丰度矩阵。")
            return
        try:
            from modules.module_metabolomics.service import MetabolomicsAnalysisService
            dataset, validation = MetabolomicsAnalysisService.inspect(source)
            self._fill_metadata(dataset.samples)
            errors = validation[validation["Status"].eq("ERROR")]
            lines = [f"适配器：{dataset.provenance.get('adapter', 'unknown')}", f"特征数：{len(dataset.abundance):,}", f"样本数：{len(dataset.samples)}", f"QC 数：{len(dataset.qc_samples)}"]
            for _, row in validation.iterrows(): lines.append(f"[{row['Status']}] {row['Check']}: {row['Message']}")
            self.log.setPlainText("\n".join(lines)); self.stage_label.setText("校验完成" if errors.empty else "存在阻断错误")
        except Exception as exc:
            QMessageBox.critical(self, "识别失败", str(exc))

    def _fill_metadata(self, frame):
        columns = ["Sample ID", "Sample type", "Group", "Factor A", "Factor B"]
        self.metadata_table.setRowCount(len(frame))
        for row_index, (_, row) in enumerate(frame.iterrows()):
            for column_index, name in enumerate(columns): self.metadata_table.setItem(row_index, column_index, QTableWidgetItem(str(row.get(name, ""))))
        self.metadata_table.resizeColumnsToContents()

    def _metadata_records(self):
        names = [self.metadata_table.horizontalHeaderItem(i).text() for i in range(self.metadata_table.columnCount())]
        records = []
        for row in range(self.metadata_table.rowCount()):
            records.append({name: (self.metadata_table.item(row, col).text().strip() if self.metadata_table.item(row, col) else "") for col, name in enumerate(names)})
        return records or None

    def _config(self):
        return {
            "input_path": self.input_edit.text().strip(), "output_parent": self.output_edit.text().strip(),
            "analysis_mode": self.mode_combo.currentData(), "factor_a": self.factor_a.text().strip() or "Group", "factor_b": self.factor_b.text().strip() or "Factor B",
            "metadata_records": self._metadata_records(), "confirmed_input_scale": self.scale_combo.currentData(),
            "missing_method": self.missing_combo.currentData(), "transformation": self.transform_combo.currentData(), "scaling": self.scaling_combo.currentData(),
            "one_factor_method": self.method_combo.currentData(), "posthoc": self.posthoc_combo.currentData(), "fdr_threshold": self.fdr_spin.value(), "fc_threshold": self.fc_spin.value(),
            "plsda_enabled": self.pls_check.isChecked(), "plsda_max_components": self.component_spin.value(), "cv_folds": self.cv_spin.value(), "permutations": self.permutation_spin.value(),
            "heatmap_top_n": self.heatmap_spin.value(), "random_seed": self.seed_spin.value(),
            "phenotype_file": self.phenotype_edit.text().strip() or None, "taxon_table": self.taxon_edit.text().strip() or None,
        }

    def start_analysis(self):
        config = self._config()
        if not config["input_path"] or not Path(config["input_path"]).exists() or not config["output_parent"]:
            QMessageBox.warning(self, "缺少路径", "请设置有效输入和输出父目录。")
            return
        if config["confirmed_input_scale"] == "unknown":
            QMessageBox.warning(self, "必须确认尺度", "正式差异分析前必须确认输入是 raw、log2 或 log10 丰度。")
            return
        Path(config["output_parent"]).mkdir(parents=True, exist_ok=True)
        handle, self.config_path = tempfile.mkstemp(prefix="sciforge_metabolomics_", suffix=".json")
        os.close(handle); Path(self.config_path).write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
        self.process = QProcess(self); self.process.setProcessChannelMode(QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read_process); self.process.finished.connect(self._finished)
        arguments = ["--run-metabolomics", "--config", self.config_path] if getattr(sys, "frozen", False) else ["-m", "modules.module_metabolomics.runner", "--config", self.config_path]
        self.process.start(sys.executable, arguments)
        self.process_buffer = ""; self.result_dir = None
        self.run_button.setEnabled(False); self.cancel_button.setEnabled(True); self.progress.setValue(0); self.log.clear(); self.stage_label.setText("启动分析进程…")

    def _read_process(self):
        if not self.process: return
        self.process_buffer += bytes(self.process.readAllStandardOutput()).decode("utf-8", errors="replace")
        lines = self.process_buffer.split("\n"); self.process_buffer = lines.pop()
        for line in lines:
            line = line.rstrip("\r")
            if line.startswith("PROGRESS\t"):
                payload = json.loads(line.split("\t", 1)[1]); self.progress.setValue(int(payload["percent"])); self.stage_label.setText(f"{payload['stage']}：{payload['message']}")
            elif line.startswith("RESULT\t"):
                payload = json.loads(line.split("\t", 1)[1]); self.result_dir = payload.get("output_dir")
                if not payload.get("success"):
                    self.log.append(f"\n错误原因：{payload.get('message', '未知错误')}\n日志文件：{payload.get('log', '')}\n建议：检查数据校验、样本分组、输入尺度和参数设置。")
                summary = payload.get("summary", {}); self.log.append("\n结果摘要：\n" + "\n".join(f"{key}: {value}" for key, value in summary.items()))
                rows = "".join(f"<tr><th>{key}</th><td>{value}</td></tr>" for key, value in summary.items())
                self.result_browsers[0].setHtml(f"<h2>Analysis Summary</h2><table>{rows}</table><p>结果目录：{self.result_dir}</p>")
            else: self.log.append(line)

    def _finished(self, exit_code, _status):
        self.run_button.setEnabled(True); self.cancel_button.setEnabled(False)
        success = exit_code == 0 and bool(self.result_dir)
        self.stage_label.setText("分析完成" if success else "分析失败或已取消")
        self.open_result_button.setEnabled(success)
        if self.config_path:
            try: Path(self.config_path).unlink(missing_ok=True)
            except OSError: pass
        self.process = None

    def cancel_analysis(self):
        if self.process and self.process.state() != QProcess.NotRunning:
            self.process.terminate()
            if not self.process.waitForFinished(3000): self.process.kill()

    def open_result(self):
        if self.result_dir: QDesktopServices.openUrl(QUrl.fromLocalFile(self.result_dir))
