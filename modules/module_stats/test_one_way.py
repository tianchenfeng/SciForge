from pathlib import Path
import sys
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common.io_utils import load_table_auto
from modules.module_stats.stats_core import run_one_way_anova


input_file = BASE_DIR / "tests" / "demo_data.xlsx"

df = load_table_auto(input_file)

result = run_one_way_anova(
    df=df,
    factor="B",
    value_col="value"
)

print("===== 单因素 ANOVA =====")
print(result["anova_table"])

print("\n===== Summary =====")
print(result["summary"])