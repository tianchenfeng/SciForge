from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
SCRIPTS_DIR = MODULE_DIR / "scripts"

TOOLBOX_VERSION = "v1.0"

WINDOW_TITLE = "16S 微生物分析软件"
MAIN_TITLE = "16S 微生物分析工具箱"
SUB_TITLE = "一键式微生物组分析流程"

SCRIPT_DISPLAY_MAP = {
    "01_qc.py": "01  质量控制",
    "02_alpha_beta.py": "02  Alpha / Beta 多样性分析",
    "03_taxonomy.py": "03  群落组成分析",
    "04_differential.py": "04  高变菌属热图",
    "05_picrust2.py": "05  PICRUSt2 功能预测",
    "06_association.py": "06  菌属-功能相关分析",
    "07_environment_correlation.py": "07  菌属-环境指标相关分析",
    "08_lefse.py": "08  候选标志菌筛选",
    "09_diff_visualization.py": "09  分层差异可视化",
    "10_export_top_genera.py": "10  主要菌群表格导出",
}

PIPELINE_ORDER = [
    "01_qc.py",
    "02_alpha_beta.py",
    "03_taxonomy.py",
    "04_differential.py",
    "05_picrust2.py",
    "06_association.py",
    "07_environment_correlation.py",
    "08_lefse.py",
    "09_diff_visualization.py",
    "10_export_top_genera.py",
]

DEFAULT_CHECKED = {
    "01_qc.py": False,
    "02_alpha_beta.py": False,
    "03_taxonomy.py": False,
    "04_differential.py": False,
    "05_picrust2.py": False,
    "06_association.py": False,
    "07_environment_correlation.py": False,
    "08_lefse.py": False,
    "09_diff_visualization.py": False,
    "10_export_top_genera.py": True,
}

FOOTER_LEFT = "Created by Tian Chengfeng"
FOOTER_RIGHT = "Email：tianchenfeng0630@gmail.com"
NOTICE_TEXT = "如遇报错，请将 _run_logs 日志发送至作者邮箱"

# ===== 新增：输入文件检查 =====
REQUIRED_INPUT_FILES = [
    "metadata.tsv",
    "feature-table.tsv",
    "taxonomy.tsv",
]

# ===== 新增：输出目录说明 =====
OUTPUT_DIR_HINTS = [
    ("_run_logs", "运行日志（报错排查用）"),

    ("01_qc", "质量控制结果"),
    ("02_alpha_beta", "Alpha / Beta 多样性分析"),
    ("03_taxa", "群落组成分析（柱状图 / 堆叠图）"),
    ("04_diff", "差异分析（高变菌）"),
    ("05_picrust2", "PICRUSt2 功能预测"),
    ("06_assoc", "菌属-功能相关分析"),
    ("07_env_corr", "菌属-环境因子相关分析"),
    ("08_candidate_biomarkers", "标志菌筛选（LEfSe）"),
    ("09_diff_viz", "分层差异可视化"),
    ("10_data_export", "主要菌群汇总表"),
]
PICRUST2_REQUIRED_FILES = [
    "KO_predicted.tsv.gz",
    "EC_predicted.tsv.gz",
]