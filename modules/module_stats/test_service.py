from pathlib import Path
import sys

BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from modules.module_stats.service import run_stats_analysis, get_table_columns
from modules.module_stats.test_batch_analysis import ensure_demo_two_way_data

input_file = ensure_demo_two_way_data()
output_dir = BASE_DIR / "output"

print("===== 表头 =====")
print(get_table_columns(input_file))

print("===== 开始分析 =====")
result = run_stats_analysis(
    input_file=input_file,
    output_dir=output_dir,
    factor1="Group",
    factor2="Time",
    value_col="pH"
)

print(result["success"])
print(result["message"])
print(result["word_path"])
print(result["excel_path"])
