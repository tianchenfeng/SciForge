from pathlib import Path
import sys

BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from common.io_utils import load_table_auto
from modules.module_stats.stats_core import run_one_way_anova
from modules.module_stats.word_export import export_one_way_word_report


input_file = BASE_DIR / "tests" / "demo_data.xlsx"
output_dir = BASE_DIR / "output"
output_dir.mkdir(exist_ok=True)

df = load_table_auto(input_file)

result = run_one_way_anova(
    df=df,
    factor="B",
    value_col="value"
)

out_path = output_dir / "one_way_stats_report.docx"

export_one_way_word_report(
    output_path=out_path,
    result=result,
    factor="B",
    value_col="value"
)

print(f"[OK] 单因素 Word 已导出: {out_path}")