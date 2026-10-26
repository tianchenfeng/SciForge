from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = BASE_DIR.parent.parent

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


from modules.module_stats.stats_core import run_one_way_anova, run_two_way_anova
from modules.module_stats.word_export import (
    export_one_way_batch_word_report,
    export_two_way_batch_word_report,
    export_word_report,
)


SUPPORTED_EXTENSIONS = [".xlsx", ".csv", ".tsv", ".txt"]


def _load_table_auto(input_file):
    input_file = Path(input_file)
    suffix = input_file.suffix.lower()

    if suffix == ".xlsx":
        return pd.read_excel(input_file)
    if suffix == ".csv":
        return pd.read_csv(input_file)
    if suffix in [".tsv", ".txt"]:
        return pd.read_csv(input_file, sep="\t")

    raise ValueError(
        f"Unsupported file extension: {suffix}. "
        f"Supported extensions: {SUPPORTED_EXTENSIONS}"
    )


def _clean_column_names(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        str(c).replace("\ufeff", "").strip()
        for c in df.columns
    ]
    if df.columns.duplicated().any():
        raise ValueError("Column names must be unique after trimming spaces / 列名去除空格后不能重复。")
    return df


def _validate_analysis_input(df, factors, metrics):
    """Reject unsafe inputs before invoking the unchanged statistical routines."""
    for factor in factors:
        labels = df[factor]
        if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
            raise ValueError(f"Factor '{factor}' has missing/blank labels / 因素列存在空白分组，请补全。")
        if labels.nunique() < 2:
            raise ValueError(f"Factor '{factor}' requires at least two levels / 因素至少需要两个水平。")
    for metric in metrics:
        values = pd.to_numeric(df[metric], errors="coerce").dropna().to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"Metric '{metric}' contains infinity / 指标含无穷大，请检查原始数据。")
        magnitudes = np.abs(values[values != 0])
        if ((magnitudes > 1e150) | (magnitudes < 1e-150)).any():
            raise ValueError(f"Metric '{metric}' exceeds the supported numeric range (1e-150 to 1e150) / 数值超出安全范围，请换算单位后重试。")


def _format_p_value(p):
    try:
        p = float(p)
    except Exception:
        return "NA"

    if pd.isna(p) or p in [float("inf"), float("-inf")]:
        return "NA"

    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def _format_number(value, digits=3):
    try:
        return f"{float(value):.{digits}f}"
    except Exception:
        return "NA"


def _format_mean_letter(row):
    try:
        value = f"{float(row['mean']):.2f}"
    except Exception:
        return "NA"

    letter = row.get("final_letter", "")
    if letter is None:
        letter = ""
    return f"{value}{letter}"


def _anova_p(anova_table, key):
    p_col = "PR(>F)"
    if p_col not in anova_table.columns:
        return None

    if key in anova_table.index:
        return anova_table.loc[key, p_col]

    key_norm = _normalize_anova_key(key)
    for index_value in anova_table.index:
        index_text = str(index_value)
        if index_text.strip().lower() == "residual":
            continue

        if _normalize_anova_key(index_text) == key_norm:
            return anova_table.loc[index_value, p_col]

    return None


def _normalize_anova_key(value):
    text = str(value).strip().replace(" ", "")
    return re.sub(r"Q\([\"'](.+?)[\"']\)", r"\1", text)


def _validate_columns(df, columns):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(missing)
            + ". Current columns: "
            + ", ".join(map(str, df.columns))
        )


def _detect_metric_columns(df, start_index):
    metrics = []

    for column in df.columns[start_index:]:
        column_name = str(column)

        numeric_values = pd.to_numeric(df[column], errors="coerce")
        if numeric_values.notna().sum() >= 2:
            metrics.append(column_name)

    return metrics


def _prepare_input(input_file, required_columns):
    input_file = Path(input_file)
    if not input_file.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_file}")

    df = _load_table_auto(input_file)
    df = _clean_column_names(df)
    _validate_columns(df, required_columns)
    return input_file, df


def _build_one_way_tables(results, factor):
    levels = []
    for result in results.values():
        for level in result["raw_df"][factor].astype(str):
            if level not in levels:
                levels.append(level)

    rows = []
    anova_rows = []
    summary_rows = []

    for metric, result in results.items():
        summary = result["summary"].copy()
        anova = result["anova_table"]
        p_factor = _anova_p(anova, f"C({factor})")

        row = {"Project": metric}
        for level in levels:
            sub = summary[summary[factor].astype(str) == str(level)]
            row[level] = "NA" if sub.empty else _format_mean_letter(sub.iloc[0])

        row["SEM"] = _format_number(summary["sem"].mean(), digits=3)
        row["P"] = _format_p_value(p_factor)
        rows.append(row)

        anova_rows.append({
            "Project": metric,
            "P": _format_p_value(p_factor),
        })

        summary = summary.copy()
        summary.insert(0, "Project", metric)
        summary_rows.append(summary)

    paper_table = pd.DataFrame(rows, columns=["Project"] + levels + ["SEM", "P"])
    anova_all = pd.DataFrame(anova_rows)
    summary_long = (
        pd.concat(summary_rows, ignore_index=True)
        if summary_rows
        else pd.DataFrame()
    )

    return paper_table, anova_all, summary_long


def _build_two_way_tables(results, factor1, factor2):
    factor1_levels = []
    factor2_levels = []

    for result in results.values():
        raw_df = result["raw_df"]
        for level in raw_df[factor1].astype(str):
            if level not in factor1_levels:
                factor1_levels.append(level)
        for level in raw_df[factor2].astype(str):
            if level not in factor2_levels:
                factor2_levels.append(level)

    rows = []
    anova_rows = []
    summary_rows = []

    for metric, result in results.items():
        summary = result["summary"].copy()
        anova = result["anova_table"]
        p_factor1 = _anova_p(anova, f"C({factor1})")
        p_factor2 = _anova_p(anova, f"C({factor2})")
        p_interaction = _anova_p(anova, f"C({factor1}):C({factor2})")

        for row_index, factor2_level in enumerate(factor2_levels):
            row = {
                "Project": metric if row_index == 0 else "",
                factor2: factor2_level,
            }

            for factor1_level in factor1_levels:
                sub = summary[
                    (summary[factor1].astype(str) == str(factor1_level))
                    & (summary[factor2].astype(str) == str(factor2_level))
                ]
                row[factor1_level] = "NA" if sub.empty else _format_mean_letter(sub.iloc[0])

            row["SEM"] = _format_number(summary["sem"].mean(), digits=3) if row_index == 0 else ""
            row[f"P_{factor1}"] = _format_p_value(p_factor1) if row_index == 0 else ""
            row[f"P_{factor2}"] = _format_p_value(p_factor2) if row_index == 0 else ""
            row[f"P_{factor1}\u00d7{factor2}"] = _format_p_value(p_interaction) if row_index == 0 else ""
            rows.append(row)

        anova_rows.append({
            "Project": metric,
            f"P_{factor1}": _format_p_value(p_factor1),
            f"P_{factor2}": _format_p_value(p_factor2),
            f"P_{factor1}\u00d7{factor2}": _format_p_value(p_interaction),
        })

        summary = summary.copy()
        summary.insert(0, "Project", metric)
        summary_rows.append(summary)

    paper_columns = (
        ["Project", factor2]
        + factor1_levels
        + ["SEM", f"P_{factor1}", f"P_{factor2}", f"P_{factor1}\u00d7{factor2}"]
    )
    paper_table = pd.DataFrame(rows, columns=paper_columns)
    anova_all = pd.DataFrame(anova_rows)
    summary_long = (
        pd.concat(summary_rows, ignore_index=True)
        if summary_rows
        else pd.DataFrame()
    )

    return paper_table, anova_all, summary_long


def _write_excel(excel_path, paper_table, anova_all, summary_long):
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        paper_table.to_excel(writer, sheet_name="paper_table", index=False)
        anova_all.to_excel(writer, sheet_name="anova_all", index=False)
        summary_long.to_excel(writer, sheet_name="summary_long", index=False)


def _build_message(base_message, paper_table):
    text = paper_table.astype(str).to_string(index=False)
    if "NA" in text:
        return base_message + " Some P values could not be estimated and are shown as NA."
    return base_message


def get_table_columns(input_file) -> list[str]:
    df = _load_table_auto(input_file)
    df = _clean_column_names(df)
    return list(df.columns)


def run_stats_analysis(input_file, output_dir, factor1, factor2, value_col):
    try:
        input_file, df = _prepare_input(input_file, [factor1, factor2, value_col])
        _validate_analysis_input(df, [factor1, factor2], [value_col])
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        result = run_two_way_anova(
            df=df,
            factor1=factor1,
            factor2=factor2,
            value_col=value_col,
        )

        base_name = input_file.stem
        word_path = output_dir / f"{base_name}_stats_report.docx"
        excel_path = output_dir / f"{base_name}_stats_summary.xlsx"

        export_word_report(
            output_path=word_path,
            result=result,
            factor1=factor1,
            factor2=factor2,
            value_col=value_col,
        )

        with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
            result["summary"].to_excel(writer, sheet_name="summary", index=False)
            result["anova_table"].to_excel(writer, sheet_name="anova")

        return {
            "success": True,
            "analysis_type": "two_way",
            "factor1": factor1,
            "factor2": factor2,
            "metrics": [value_col],
            "results": {value_col: result},
            "paper_table": None,
            "word_path": str(word_path),
            "excel_path": str(excel_path),
            "message": "Statistics analysis completed.",
            "summary": result["summary"],
            "anova_table": result["anova_table"],
        }

    except Exception as e:
        return {
            "success": False,
            "analysis_type": "two_way",
            "factor1": factor1,
            "factor2": factor2,
            "metrics": [],
            "results": {},
            "paper_table": None,
            "word_path": "",
            "excel_path": "",
            "message": str(e),
        }


def run_one_way_batch_analysis(input_file, output_dir, factor=None):
    try:
        input_file, df = _prepare_input(input_file, [])
        if len(df.columns) < 2:
            raise ValueError("One-way input requires one factor column and at least one metric column.")

        factor = str(df.columns[0])
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        metrics = _detect_metric_columns(df, start_index=1)
        if not metrics:
            raise ValueError("No metric columns found with at least two numeric values.")

        _validate_analysis_input(df, [factor], metrics)

        results = {}
        for metric in metrics:
            results[metric] = run_one_way_anova(
                df=df,
                factor=factor,
                value_col=metric,
            )

        paper_table, anova_all, summary_long = _build_one_way_tables(results, factor)

        base_name = input_file.stem
        word_path = output_dir / f"{base_name}_one_way_batch_report.docx"
        excel_path = output_dir / f"{base_name}_one_way_batch_summary.xlsx"

        batch_result = {
            "success": True,
            "analysis_type": "one_way",
            "factor1": factor,
            "factor2": None,
            "metrics": metrics,
            "results": results,
            "paper_table": paper_table,
            "anova_all": anova_all,
            "summary_long": summary_long,
            "word_path": str(word_path),
            "excel_path": str(excel_path),
            "message": _build_message("One-way batch analysis completed.", paper_table),
        }

        export_one_way_batch_word_report(word_path, batch_result, factor)
        _write_excel(excel_path, paper_table, anova_all, summary_long)

        return batch_result

    except Exception as e:
        return {
            "success": False,
            "analysis_type": "one_way",
            "factor1": factor,
            "factor2": None,
            "metrics": [],
            "results": {},
            "paper_table": None,
            "anova_all": None,
            "summary_long": None,
            "word_path": "",
            "excel_path": "",
            "message": str(e),
        }


def run_two_way_batch_analysis(input_file, output_dir, factor1=None, factor2=None):
    try:
        input_file, df = _prepare_input(input_file, [])
        if len(df.columns) < 3:
            raise ValueError("Two-way input requires two factor columns and at least one metric column.")

        factor1 = str(df.columns[0])
        factor2 = str(df.columns[1])
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        metrics = _detect_metric_columns(df, start_index=2)
        if not metrics:
            raise ValueError("No metric columns found with at least two numeric values.")

        _validate_analysis_input(df, [factor1, factor2], metrics)

        results = {}
        for metric in metrics:
            results[metric] = run_two_way_anova(
                df=df,
                factor1=factor1,
                factor2=factor2,
                value_col=metric,
            )

        paper_table, anova_all, summary_long = _build_two_way_tables(results, factor1, factor2)

        base_name = input_file.stem
        word_path = output_dir / f"{base_name}_two_way_batch_report.docx"
        excel_path = output_dir / f"{base_name}_two_way_batch_summary.xlsx"

        batch_result = {
            "success": True,
            "analysis_type": "two_way",
            "factor1": factor1,
            "factor2": factor2,
            "metrics": metrics,
            "results": results,
            "paper_table": paper_table,
            "anova_all": anova_all,
            "summary_long": summary_long,
            "word_path": str(word_path),
            "excel_path": str(excel_path),
            "message": _build_message("Two-way batch analysis completed.", paper_table),
        }

        export_two_way_batch_word_report(word_path, batch_result, factor1, factor2)
        _write_excel(excel_path, paper_table, anova_all, summary_long)

        return batch_result

    except Exception as e:
        return {
            "success": False,
            "analysis_type": "two_way",
            "factor1": factor1,
            "factor2": factor2,
            "metrics": [],
            "results": {},
            "paper_table": None,
            "anova_all": None,
            "summary_long": None,
            "word_path": "",
            "excel_path": "",
            "message": str(e),
        }
