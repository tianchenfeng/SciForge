from __future__ import annotations

from pathlib import Path

import pandas as pd

from modules.module_gas_fit.export_results import export_excel
from modules.module_gas_fit.fitting import fit_sample
from modules.module_gas_fit.models import MODEL_EXPONENTIAL, MODEL_LAG
from modules.module_gas_fit.plotting import draw_fit_plot, draw_summary_plot
from modules.module_gas_fit.preprocess import clean_gas_data, load_gas_table


class GasFitService:
    @staticmethod
    def run(input_file: str, output_dir: str, model_name: str = MODEL_EXPONENTIAL, logger=None) -> dict:
        if model_name not in {MODEL_EXPONENTIAL, MODEL_LAG}:
            raise ValueError("model_name must be exponential or lag")

        input_path = Path(input_file)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        def log(message: str):
            if callable(logger):
                logger(message)

        log(f"Reading input: {input_path}")
        df = load_gas_table(input_path)
        cleaned = clean_gas_data(df)
        data = cleaned["df"]
        sample_col = cleaned["sample_col"]
        time_columns = cleaned["time_columns"]
        times = cleaned["times"]
        log(f"Detected time columns: {times}")

        success_results = []
        result_rows = []
        failed_rows = []
        for _, row in data.iterrows():
            sample = str(row[sample_col])
            values = pd.to_numeric(row[time_columns], errors="coerce").to_numpy(dtype=float)
            try:
                result = fit_sample(sample, times, values, model_name)
            except Exception as exc:
                error_message = str(exc)
                failed_row = {
                    "Sample": sample,
                    "Model": model_name,
                    "A": None,
                    "c": None,
                    "Lag": None,
                    "R2": None,
                    "RMSE": None,
                    "status": "failed",
                    "error_message": error_message,
                }
                result_rows.append(failed_row)
                failed_rows.append(failed_row)
                log(f"Failed {sample}: {error_message}")
                continue

            success_results.append(result)
            row_result = result.to_row()
            row_result["status"] = "success"
            row_result["error_message"] = ""
            result_rows.append(row_result)
            log(f"Fitted {sample}: A={result.a:.3f}, c={result.c:.4f}, Lag={result.lag:.3f}, R2={result.r2:.4f}")

        plot_dir = output_path / "plots"
        plot_files = [draw_fit_plot(result, plot_dir / f"{_safe_name(result.sample)}.png") for result in success_results]
        summary_plot = draw_summary_plot(success_results, output_path / "gas_fit_summary.png") if success_results else None
        excel_file = export_excel(result_rows, failed_rows, success_results, output_path / "gas_fit_results.xlsx")
        log("Gas fitting finished")

        return {
            "success": True,
            "message": "Gas fitting finished",
            "excel": excel_file,
            "plots": plot_files,
            "summary_plot": summary_plot,
            "results": result_rows,
            "failed_samples": failed_rows,
            "success_count": len(success_results),
            "failed_count": len(failed_rows),
        }


def _safe_name(value: str) -> str:
    keep = []
    for ch in str(value):
        keep.append(ch if ch.isalnum() or ch in ("-", "_") else "_")
    return "".join(keep).strip("_") or "sample"
