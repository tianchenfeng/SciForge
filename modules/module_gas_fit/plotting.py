from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from modules.module_gas_fit.models import predict


def _font(size: int = 14):
    try:
        return ImageFont.truetype("arial.ttf", size)
    except OSError:
        return ImageFont.load_default()


def _scale(points, x_min, x_max, y_min, y_max, box):
    left, top, right, bottom = box
    scaled = []
    for x, y in points:
        sx = left + (x - x_min) / (x_max - x_min or 1.0) * (right - left)
        sy = bottom - (y - y_min) / (y_max - y_min or 1.0) * (bottom - top)
        scaled.append((float(sx), float(sy)))
    return scaled


def draw_fit_plot(result, output_path: str | Path) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    width, height = 900, 560
    box = (80, 60, 850, 470)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    t = np.asarray(result.times, dtype=float)
    y = np.asarray(result.observed, dtype=float)
    x_min, x_max = float(np.min(t)), float(np.max(t))
    y_min = min(0.0, float(np.min(y)))
    y_max = max(float(np.max(y)), float(np.max(result.fitted))) * 1.12 or 1.0

    draw.rectangle(box, outline="#334155", width=2)
    draw.text((80, 24), f"{result.sample} gas fit ({result.model})", fill="#111827", font=_font(18))
    draw.text((80, 500), "Time", fill="#111827", font=_font())
    draw.text((12, 245), "Gas", fill="#111827", font=_font())

    for i in range(6):
        x = box[0] + i * (box[2] - box[0]) / 5
        y_grid = box[3] - i * (box[3] - box[1]) / 5
        draw.line((x, box[1], x, box[3]), fill="#e5e7eb")
        draw.line((box[0], y_grid, box[2], y_grid), fill="#e5e7eb")

    obs_points = _scale(zip(t, y), x_min, x_max, y_min, y_max, box)
    dense_t = np.linspace(x_min, x_max, 240)
    dense_y = predict(result.model, dense_t, result.a, result.c, result.lag)
    fit_points = _scale(zip(dense_t, dense_y), x_min, x_max, y_min, y_max, box)

    if len(fit_points) >= 2:
        draw.line(fit_points, fill="#2563eb", width=3)
    for x, yv in obs_points:
        draw.ellipse((x - 4, yv - 4, x + 4, yv + 4), fill="#dc2626", outline="#991b1b")

    legend = f"A={result.a:.3f}, c={result.c:.4f}, Lag={result.lag:.3f}, R2={result.r2:.4f}, RMSE={result.rmse:.3f}"
    draw.text((330, 500), legend, fill="#334155", font=_font())
    image.save(output_path)
    return str(output_path)


def draw_summary_plot(results, output_path: str | Path) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 1000, 620
    box = (80, 60, 940, 520)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    all_t = [x for result in results for x in result.times]
    all_y = [x for result in results for x in result.observed + result.fitted]
    x_min, x_max = min(all_t), max(all_t)
    y_min, y_max = min(0.0, min(all_y)), max(all_y) * 1.12 or 1.0
    colors = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#ea580c", "#0891b2", "#4f46e5"]

    draw.rectangle(box, outline="#334155", width=2)
    draw.text((80, 24), "Gas fit summary", fill="#111827", font=_font(18))
    for idx, result in enumerate(results):
        dense_t = np.linspace(min(result.times), max(result.times), 200)
        dense_y = predict(result.model, dense_t, result.a, result.c, result.lag)
        points = _scale(zip(dense_t, dense_y), x_min, x_max, y_min, y_max, box)
        draw.line(points, fill=colors[idx % len(colors)], width=2)
        if idx < 12:
            draw.text((80 + (idx % 4) * 220, 540 + (idx // 4) * 20), result.sample, fill=colors[idx % len(colors)], font=_font())

    image.save(output_path)
    return str(output_path)
