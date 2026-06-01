from pathlib import Path
import sys
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common.io_utils import load_table_auto, find_first_table
from stats_core import run_two_way_anova
from word_export import export_word_report


TEST_DIR = BASE_DIR / "tests"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)


REQUIRED_COLUMNS = ["Sample", "A", "B", "value"]


def check_input_columns(df: pd.DataFrame):
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]

    if missing:
        raise ValueError(
            "输入表缺少必要列："
            + ", ".join(missing)
            + "\n当前表格列名为："
            + ", ".join(map(str, df.columns))
        )


def main():
    print("===== 读取数据 =====")

    input_file = find_first_table(TEST_DIR)
    print(f"检测到输入文件: {input_file}")

    df = load_table_auto(input_file)
    check_input_columns(df)

    print(df.head())

    print("\n===== 数据摘要 =====")
    print(f"样本数: {df.shape[0]}")
    print(f"A 水平: {', '.join(map(str, sorted(df['A'].unique())))}")
    print(f"B 水平: {', '.join(map(str, sorted(df['B'].astype(str).unique())))}")

    print("\n===== 双因素分析 =====")

    result = run_two_way_anova(
        df=df,
        factor1="A",
        factor2="B",
        value_col="value"
    )

    print(result["anova_table"])

    word_path = OUTPUT_DIR / "stats_report.docx"
    excel_path = OUTPUT_DIR / "stats_summary.xlsx"

    export_word_report(
        output_path=word_path,
        result=result,
        factor1="A",
        factor2="B",
        value_col="value"
    )

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        result["summary"].to_excel(writer, sheet_name="summary", index=False)
        result["anova_table"].to_excel(writer, sheet_name="anova")

    print("\n===== 输出完成 =====")
    print(f"Word:  {word_path}")
    print(f"Excel: {excel_path}")


if __name__ == "__main__":
    main()