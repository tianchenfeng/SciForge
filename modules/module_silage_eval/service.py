from __future__ import annotations

from pathlib import Path

import pandas as pd

from modules.module_silage_eval.entropy_weight import entropy_weights
from modules.module_silage_eval.membership import infer_indicator_directions, membership_scores, range_normalize
from modules.module_silage_eval.ranking import build_ranking
from modules.module_silage_eval.word_export import export_word_report


class SilageEvalService:
    @staticmethod
    def inspect_directions(input_file: str) -> dict[str, str]:
        df = _read_table(input_file)
        if df.shape[1] < 2:
            raise ValueError("Input table must contain Treatment and indicator columns")
        _, indicator_cols = _select_numeric_indicators(df)
        return infer_indicator_directions(indicator_cols)

    @staticmethod
    def run(input_file: str, output_dir: str, directions: dict[str, str] | None = None, logger=None) -> dict:
        input_path = Path(input_file)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        def log(message: str):
            if callable(logger):
                logger(message)

        log(f"Reading input: {input_path}")
        df = _read_table(input_path)
        if df.shape[1] < 2:
            raise ValueError("Input table must contain Treatment and indicator columns")

        treatment_col = df.columns[0]
        numeric, indicator_cols = _select_numeric_indicators(df)
        if not indicator_cols:
            raise ValueError("No numeric indicator data found")

        treatment_order = df[treatment_col].dropna().astype(str).drop_duplicates().tolist()
        treatment_mean = pd.concat([df[[treatment_col]].astype(str), numeric[indicator_cols]], axis=1)
        treatment_mean = treatment_mean.groupby(treatment_col, sort=False, as_index=False).mean(numeric_only=True)
        treatment_mean[treatment_col] = pd.Categorical(treatment_mean[treatment_col], categories=treatment_order, ordered=True)
        treatment_mean = treatment_mean.sort_values(treatment_col).reset_index(drop=True)
        treatment_mean[treatment_col] = treatment_mean[treatment_col].astype(str)

        direction_map = directions or infer_indicator_directions(indicator_cols)
        direction_map = {col: direction_map.get(col, infer_indicator_directions([col])[col]) for col in indicator_cols}

        log(f"Detected {len(indicator_cols)} numeric indicators")
        normalized = range_normalize(treatment_mean[indicator_cols], direction_map)
        weights = entropy_weights(normalized)
        scores = membership_scores(normalized, weights)
        ranking = build_ranking(treatment_mean[treatment_col], scores)
        weight_df = pd.DataFrame({"Indicator": weights.index, "Weight": weights.values})
        direction_df = pd.DataFrame({"Indicator": list(direction_map.keys()), "Direction": list(direction_map.values())})
        normalized_export = pd.concat([treatment_mean[[treatment_col]], normalized], axis=1)

        excel_file = output_path / "silage_eval_results.xlsx"
        with pd.ExcelWriter(excel_file, engine="openpyxl") as writer:
            df.to_excel(writer, index=False, sheet_name="raw_data")
            treatment_mean.to_excel(writer, index=False, sheet_name="treatment_mean")
            direction_df.to_excel(writer, index=False, sheet_name="directions")
            normalized_export.to_excel(writer, index=False, sheet_name="normalized")
            weight_df.to_excel(writer, index=False, sheet_name="entropy_weights")
            ranking.to_excel(writer, index=False, sheet_name="ranking")

        word_file = export_word_report(ranking, weight_df, direction_map, output_path / "silage_eval_report.docx")
        log("Silage evaluation finished")
        return {
            "success": True,
            "message": "Silage evaluation finished",
            "excel": str(excel_file),
            "word": word_file,
            "directions": direction_map,
            "treatment_mean": treatment_mean.to_dict(orient="records"),
            "ranking": ranking.to_dict(orient="records"),
        }


def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    return pd.read_csv(path)


def _select_numeric_indicators(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    numeric = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    indicator_cols = [col for col in numeric.columns if numeric[col].notna().sum() >= 2]
    return numeric, indicator_cols
