from __future__ import annotations

from pathlib import Path

import pandas as pd


COLUMN_ALIASES = {
    "gene": ["Gene", "gene", "gene_id", "GeneID", "Symbol", "symbol", "id"],
    "log2fc": ["log2FC", "log2FoldChange", "logFC", "Log2FC"],
    "pvalue": ["pvalue", "PValue", "p_value", "P.Value", "P"],
    "padj": ["padj", "FDR", "qvalue", "q_value", "adj.P.Val", "qvalue"],
}


def read_table(input_file: str | Path) -> pd.DataFrame:
    path = Path(input_file)
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")
    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def read_deg_table(input_file: str | Path) -> pd.DataFrame:
    return read_table(input_file)


def detect_columns(df: pd.DataFrame) -> dict[str, str]:
    detected = {}
    lower_map = {str(col).lower(): col for col in df.columns}
    for key, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias.lower() in lower_map:
                detected[key] = lower_map[alias.lower()]
                break
    missing = [key for key in COLUMN_ALIASES if key not in detected]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    return detected


def detect_gene_column(df: pd.DataFrame) -> str:
    lower_map = {str(col).lower(): col for col in df.columns}
    for alias in COLUMN_ALIASES["gene"]:
        if alias.lower() in lower_map:
            return lower_map[alias.lower()]
    return df.columns[0]


def standardize_deg_table(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    cols = detect_columns(df)
    result = pd.DataFrame({
        "Gene": df[cols["gene"]].astype(str),
        "log2FC": pd.to_numeric(df[cols["log2fc"]], errors="coerce"),
        "pvalue": pd.to_numeric(df[cols["pvalue"]], errors="coerce"),
        "padj": pd.to_numeric(df[cols["padj"]], errors="coerce"),
    })
    for col in df.columns:
        if col not in cols.values() and col not in result.columns:
            result[col] = df[col]
    result = result.dropna(subset=["Gene", "log2FC", "padj"]).reset_index(drop=True)
    return result, cols


def load_annotation(annotation_file: str | Path | None) -> pd.DataFrame:
    if not annotation_file:
        return pd.DataFrame(columns=["Gene"])
    annotation = read_table(annotation_file)
    gene_col = detect_gene_column(annotation)
    annotation = annotation.rename(columns={gene_col: "Gene"}).copy()
    annotation["Gene"] = annotation["Gene"].astype(str)
    return annotation.drop_duplicates(subset=["Gene"])


def load_expression_matrix(expression_file: str | Path | None) -> pd.DataFrame:
    if not expression_file:
        return pd.DataFrame()
    expression = read_table(expression_file)
    gene_col = detect_gene_column(expression)
    expression = expression.rename(columns={gene_col: "Gene"}).copy()
    expression["Gene"] = expression["Gene"].astype(str)
    for col in expression.columns:
        if col != "Gene":
            expression[col] = pd.to_numeric(expression[col], errors="coerce")
    return expression.dropna(subset=["Gene"]).drop_duplicates(subset=["Gene"])
