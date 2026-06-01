from pathlib import Path
from service import run_stats_analysis, get_table_columns

BASE_DIR = Path(__file__).resolve().parent

input_file = BASE_DIR / "tests" / "demo_data.xlsx"
output_dir = BASE_DIR / "output"

print("===== 表头 =====")
print(get_table_columns(input_file))

print("===== 开始分析 =====")
result = run_stats_analysis(
    input_file=input_file,
    output_dir=output_dir,
    factor1="A",
    factor2="B",
    value_col="value"
)

print(result["success"])
print(result["message"])
print(result["word_path"])
print(result["excel_path"])