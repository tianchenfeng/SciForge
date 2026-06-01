from __future__ import annotations

from pathlib import Path

import pandas as pd


FIT_COLUMNS = ["Sample", "Model", "A", "c", "Lag", "R2", "RMSE", "status", "error_message"]
FAILED_COLUMNS = FIT_COLUMNS
FITTED_VALUE_COLUMNS = ["Sample", "Time", "Observed", "Fitted", "Residual"]


def export_excel(result_rows, failed_rows, success_results, output_path: str | Path) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fit_df = pd.DataFrame(result_rows, columns=FIT_COLUMNS)
    failed_df = pd.DataFrame(failed_rows, columns=FAILED_COLUMNS)
    long_rows = []
    for result in success_results:
        for t, obs, fit in zip(result.times, result.observed, result.fitted):
            long_rows.append({
                "Sample": result.sample,
                "Time": t,
                "Observed": obs,
                "Fitted": fit,
                "Residual": obs - fit,
            })
    long_df = pd.DataFrame(long_rows, columns=FITTED_VALUE_COLUMNS)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        fit_df.to_excel(writer, index=False, sheet_name="fit_results")
        failed_df.to_excel(writer, index=False, sheet_name="failed_samples")
        long_df.to_excel(writer, index=False, sheet_name="fitted_values")
    return str(output_path)
