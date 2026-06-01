from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_gas_table(input_file: str | Path) -> pd.DataFrame:
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    suffix = input_path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        df = pd.read_excel(input_path)
    elif suffix in {".csv", ".txt"}:
        df = pd.read_csv(input_path)
    else:
        raise ValueError(f"Unsupported input file type: {suffix}")

    if df.empty:
        raise ValueError("Input table is empty")
    if df.shape[1] < 2:
        raise ValueError("Input table must contain a sample column and time columns")
    return df


def parse_time_columns(df: pd.DataFrame) -> list[tuple[str, float]]:
    time_columns = []
    for col in df.columns[1:]:
        try:
            time_value = float(str(col).strip())
        except (TypeError, ValueError):
            continue
        time_columns.append((col, time_value))

    if len(time_columns) < 3:
        raise ValueError("At least 3 parseable time columns are required")
    return time_columns


def clean_gas_data(df: pd.DataFrame) -> dict:
    sample_col = df.columns[0]
    time_pairs = parse_time_columns(df)
    time_columns = [col for col, _ in time_pairs]
    times = [time_value for _, time_value in time_pairs]

    cleaned = df[[sample_col] + time_columns].copy()
    cleaned[sample_col] = cleaned[sample_col].fillna("").astype(str).str.strip()
    empty_sample_mask = cleaned[sample_col] == ""
    for seq, idx in enumerate(cleaned.index[empty_sample_mask], start=1):
        cleaned.at[idx, sample_col] = f"sample_{seq}"
    for col in time_columns:
        cleaned[col] = pd.to_numeric(cleaned[col], errors="coerce")

    return {
        "df": cleaned,
        "sample_col": sample_col,
        "time_columns": time_columns,
        "times": times,
    }
