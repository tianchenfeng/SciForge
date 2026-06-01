from __future__ import annotations

import numpy as np


MODEL_EXPONENTIAL = "exponential"
MODEL_LAG = "lag"


def exponential_model(t, a: float, c: float):
    t_arr = np.asarray(t, dtype=float)
    return a * (1.0 - np.exp(-c * t_arr))


def lag_model(t, a: float, c: float, lag: float):
    t_arr = np.asarray(t, dtype=float)
    effective_t = np.maximum(t_arr - lag, 0.0)
    return a * (1.0 - np.exp(-c * effective_t))


def predict(model_name: str, t, a: float, c: float, lag: float = 0.0):
    if model_name == MODEL_EXPONENTIAL:
        return exponential_model(t, a, c)
    if model_name == MODEL_LAG:
        return lag_model(t, a, c, lag)
    raise ValueError(f"Unsupported gas model: {model_name}")
