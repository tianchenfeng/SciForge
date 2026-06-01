from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def draw_volcano_plot(df: pd.DataFrame, output_path: str | Path, log2fc_threshold: float, padj_threshold: float) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_df = df.copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
    colors = np.where(
        (plot_df["log2FC"].abs() >= log2fc_threshold) & (plot_df["padj"] <= padj_threshold),
        np.where(plot_df["log2FC"] > 0, "#dc2626", "#2563eb"),
        "#94a3b8",
    )

    plt.figure(figsize=(8.5, 6))
    plt.scatter(plot_df["log2FC"], plot_df["neg_log10_padj"], c=colors, s=24, alpha=0.85, edgecolors="none")
    plt.axvline(log2fc_threshold, color="#475569", linestyle="--", linewidth=1)
    plt.axvline(-log2fc_threshold, color="#475569", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(padj_threshold), color="#475569", linestyle="--", linewidth=1)
    plt.xlabel("log2FC")
    plt.ylabel("-log10(padj)")
    plt.title("Volcano plot")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    return str(output_path)


def draw_enrichment_barplot(df: pd.DataFrame, output_path: str | Path, title: str, top_n: int = 10) -> str:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plot_df = df.head(top_n).copy()
    plt.figure(figsize=(9, 5.5))
    if plot_df.empty:
        plt.text(0.5, 0.5, "No enrichment results", ha="center", va="center")
        plt.axis("off")
    else:
        plot_df["score"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
        labels = plot_df["Term"].astype(str).str.slice(0, 60)
        plt.barh(labels[::-1], plot_df["score"][::-1], color="#4f46e5")
        plt.xlabel("-log10(padj)")
        plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    return str(output_path)


def draw_candidate_heatmap(expression_df: pd.DataFrame, candidate_df: pd.DataFrame, output_path: str | Path, top_n: int = 20) -> str | None:
    if expression_df.empty or candidate_df.empty:
        return None
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    genes = candidate_df["Gene"].astype(str).head(top_n).tolist()
    matrix = expression_df[expression_df["Gene"].astype(str).isin(genes)].copy()
    if matrix.empty:
        return None
    matrix = matrix.set_index("Gene")
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    matrix = matrix.dropna(how="all")
    if matrix.empty:
        return None
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, 1), axis=0).fillna(0)

    plt.figure(figsize=(max(6, z.shape[1] * 0.8), max(4, z.shape[0] * 0.35)))
    plt.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    plt.colorbar(label="row z-score")
    plt.yticks(range(len(z.index)), z.index)
    plt.xticks(range(len(z.columns)), z.columns, rotation=45, ha="right")
    plt.title("Candidate gene heatmap")
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    return str(output_path)
