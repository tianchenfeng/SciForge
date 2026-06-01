from pathlib import Path
import sys

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from modules.module_stats.service import (
    run_one_way_batch_analysis,
    run_two_way_batch_analysis,
)


def ensure_demo_one_way_data():
    tests_dir = BASE_DIR / "tests"
    tests_dir.mkdir(exist_ok=True)
    input_file = tests_dir / "demo_one_way.xlsx"

    rows = []
    for b_index, b_level in enumerate([3, 7, 15, 30]):
        for rep in range(1, 4):
            base = 10 + b_index * 1.7 + rep * 0.25
            rows.append({
                "B": b_level,
                "value": base,
                "DML": 1.8 + b_index * 0.35 + rep * 0.05,
                "PH": 4.8 - b_index * 0.12 + rep * 0.01,
                "Remark": "text note",
                "EmptyCol": None,
            })

    pd.DataFrame(rows).to_excel(input_file, index=False)
    return input_file


def ensure_demo_two_way_data():
    tests_dir = BASE_DIR / "tests"
    tests_dir.mkdir(exist_ok=True)
    input_file = tests_dir / "demo_two_way.xlsx"

    rows = []
    groups = ["C", "L", "M"]
    times = [3, 7, 15, 30]

    for group_index, group in enumerate(groups):
        for time_index, time in enumerate(times):
            for rep in range(1, 4):
                base = 10 + group_index * 2 + time_index * 1.5 + rep * 0.2
                rows.append({
                    "Group": group,
                    "Time": time,
                    "pH": 6.5 - base * 0.03,
                    "NH3-N": base + 0.5,
                    "CP%": 12 + group_index * 1.1 + rep * 0.1,
                    "LA (g/kg)": base * 1.8,
                    "DML%": 40 + group_index * 2 + time_index * 1.2,
                    "Remark": "text note",
                    "EmptyCol": None,
                })

    pd.DataFrame(rows).to_excel(input_file, index=False)
    return input_file


def _assert_output(result):
    assert result["success"], result["message"]
    assert result["metrics"], "No metrics were detected."
    assert result["paper_table"] is not None, "paper_table is missing."
    assert result["anova_all"] is not None, "anova_all is missing."
    assert result["summary_long"] is not None, "summary_long is missing."
    assert Path(result["word_path"]).exists(), result["word_path"]
    assert Path(result["excel_path"]).exists(), result["excel_path"]


def _assert_no_unrelated_or_nan(result):
    paper_text = result["paper_table"].astype(str).to_string(index=False).lower()
    anova_text = result["anova_all"].astype(str).to_string(index=False).lower()
    paper_columns = list(result["paper_table"].columns)

    assert "Remark" not in paper_columns
    assert "EmptyCol" not in paper_columns
    assert "remark" not in paper_text
    assert "emptycol" not in paper_text
    assert "nan" not in paper_text
    assert "nan" not in anova_text


def test_one_way_batch_analysis():
    input_file = ensure_demo_one_way_data()
    output_dir = BASE_DIR / "output" / "batch_test"

    result = run_one_way_batch_analysis(
        input_file=input_file,
        output_dir=output_dir,
        factor="ignored_by_fixed_format",
    )

    _assert_output(result)
    _assert_no_unrelated_or_nan(result)
    assert result["factor1"] == "B"
    assert result["metrics"] == ["value", "DML", "PH"]
    assert "B" not in result["metrics"]
    assert "Remark" not in result["metrics"]
    assert "EmptyCol" not in result["metrics"]
    assert list(result["paper_table"].columns[:5]) == ["Project", "3", "7", "15", "30"]

    print("one_way_metrics=", result["metrics"])
    print("one_way_group_order=", list(result["paper_table"].columns[1:-2]))
    print("one_way_word_path=", result["word_path"])
    print("one_way_excel_path=", result["excel_path"])


def test_two_way_batch_analysis():
    input_file = ensure_demo_two_way_data()
    output_dir = BASE_DIR / "output" / "batch_test"

    result = run_two_way_batch_analysis(
        input_file=input_file,
        output_dir=output_dir,
        factor1="ignored_a",
        factor2="ignored_b",
    )

    _assert_output(result)
    _assert_no_unrelated_or_nan(result)
    assert result["factor1"] == "Group"
    assert result["factor2"] == "Time"
    assert result["metrics"] == ["pH", "NH3-N", "CP%", "LA (g/kg)", "DML%"]
    assert "Time" not in result["metrics"]
    assert "Remark" not in result["metrics"]
    assert "EmptyCol" not in result["metrics"]
    assert list(result["paper_table"].columns[2:5]) == ["C", "L", "M"]
    assert list(result["paper_table"]["Time"].drop_duplicates().astype(str)) == ["3", "7", "15", "30"]

    print("two_way_metrics=", result["metrics"])
    print("two_way_factor1_order=", list(result["paper_table"].columns[2:-4]))
    print("two_way_factor2_order=", list(result["paper_table"]["Time"].drop_duplicates().astype(str)))
    print("two_way_word_path=", result["word_path"])
    print("two_way_excel_path=", result["excel_path"])


if __name__ == "__main__":
    test_one_way_batch_analysis()
    test_two_way_batch_analysis()
    print("[OK] Batch analysis tests passed.")
