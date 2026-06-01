from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from modules.module_gas_fit.models import MODEL_EXPONENTIAL, MODEL_LAG, predict


@dataclass
class FitResult:
    sample: str
    model: str
    a: float
    c: float
    lag: float
    r2: float
    rmse: float
    times: list[float]
    observed: list[float]
    fitted: list[float]

    def to_row(self) -> dict:
        return {
            "Sample": self.sample,
            "Model": self.model,
            "A": self.a,
            "c": self.c,
            "Lag": self.lag,
            "R2": self.r2,
            "RMSE": self.rmse,
        }


def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    residual = y_true - y_pred
    ss_res = float(np.sum(residual ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    return r2, rmse


def _fit_for_fixed_shape(y: np.ndarray, shape: np.ndarray) -> tuple[float, float]:
    denom = float(np.sum(shape ** 2))
    if denom <= 0:
        return 0.0, float("inf")
    a = max(float(np.sum(y * shape) / denom), 0.0)
    sse = float(np.sum((y - a * shape) ** 2))
    return a, sse


def _candidate_rates(times: np.ndarray) -> np.ndarray:
    positive = times[times > 0]
    scale = float(np.max(positive)) if positive.size else 48.0
    lower = max(0.001, 0.05 / scale)
    upper = max(0.01, 8.0 / scale)
    return np.unique(np.concatenate([
        np.linspace(lower, upper, 120),
        np.geomspace(lower, upper, 80),
    ]))


def fit_sample(sample: str, times, values, model_name: str) -> FitResult:
    t = np.asarray(times, dtype=float)
    y = np.asarray(values, dtype=float)
    valid = np.isfinite(t) & np.isfinite(y)
    t = t[valid]
    y = y[valid]
    if t.size < 3:
        raise ValueError(f"{sample}: at least 3 valid time points are required")

    rates = _candidate_rates(t)
    best = {"sse": float("inf"), "a": 0.0, "c": 0.0, "lag": 0.0}

    if model_name == MODEL_EXPONENTIAL:
        for c in rates:
            shape = 1.0 - np.exp(-c * t)
            a, sse = _fit_for_fixed_shape(y, shape)
            if sse < best["sse"]:
                best.update({"sse": sse, "a": a, "c": float(c), "lag": 0.0})
    elif model_name == MODEL_LAG:
        max_lag = max(0.0, float(np.max(t)) * 0.5)
        lag_grid = np.linspace(0.0, max_lag, 61)
        for lag in lag_grid:
            effective_t = np.maximum(t - lag, 0.0)
            for c in rates:
                shape = 1.0 - np.exp(-c * effective_t)
                a, sse = _fit_for_fixed_shape(y, shape)
                if sse < best["sse"]:
                    best.update({"sse": sse, "a": a, "c": float(c), "lag": float(lag)})
    else:
        raise ValueError(f"Unsupported gas model: {model_name}")

    fitted = predict(model_name, t, best["a"], best["c"], best["lag"])
    r2, rmse = _metrics(y, fitted)
    return FitResult(
        sample=str(sample),
        model=model_name,
        a=float(best["a"]),
        c=float(best["c"]),
        lag=float(best["lag"]),
        r2=float(r2),
        rmse=float(rmse),
        times=t.tolist(),
        observed=y.tolist(),
        fitted=fitted.tolist(),
    )
