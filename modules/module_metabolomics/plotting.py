from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage


plt.rcParams.update({"font.family": ["Microsoft YaHei", "DejaVu Sans", "Arial"], "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})
GROUP_COLORS = ["#4f46e5", "#dc2626", "#059669", "#d97706", "#0891b2", "#7c3aed", "#475569"]


def save_figure(fig, base_path: str | Path, dpi: int = 600) -> dict[str, str]:
    base = Path(base_path); base.parent.mkdir(parents=True, exist_ok=True)
    paths = {}
    for extension in ["png", "pdf", "svg"]:
        path = base.with_suffix(f".{extension}")
        fig.savefig(path, dpi=dpi if extension == "png" else None, bbox_inches="tight")
        paths[extension] = str(path)
    plt.close(fig)
    return paths


def _ellipse(ax, x, y, color):
    if len(x) < 3: return
    covariance = np.cov(x, y)
    if not np.isfinite(covariance).all() or np.linalg.matrix_rank(covariance) < 2: return
    values, vectors = np.linalg.eigh(covariance); order = values.argsort()[::-1]
    values, vectors = values[order], vectors[:, order]
    angle = np.degrees(np.arctan2(*vectors[:, 0][::-1]))
    width, height = 2 * np.sqrt(values * 5.991)
    ax.add_patch(Ellipse((np.mean(x), np.mean(y)), width, height, angle=angle, fill=False, color=color, lw=1.2, alpha=0.8))


def plot_missing(missing_sample, missing_feature, base):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].bar(np.arange(len(missing_sample)), missing_sample["Missing fraction"], color="#4f46e5")
    axes[0].set_xticks(np.arange(len(missing_sample))); axes[0].set_xticklabels(missing_sample["Sample ID"], rotation=90)
    axes[0].set_ylabel("Missing fraction"); axes[0].set_title("Missingness by sample")
    axes[1].hist(missing_feature["Missing fraction"].dropna(), bins=30, color="#7c3aed", edgecolor="white")
    axes[1].set_xlabel("Missing fraction"); axes[1].set_ylabel("Features"); axes[1].set_title("Feature missingness distribution")
    fig.tight_layout(); return save_figure(fig, base)


def plot_intensity(before, after, base):
    fig, axes = plt.subplots(2, 1, figsize=(max(9, before.shape[1] * .35), 7), sharex=True)
    for ax, matrix, title in [(axes[0], before, "Before preprocessing"), (axes[1], after, "After preprocessing")]:
        ax.boxplot([pd.to_numeric(matrix[c], errors="coerce").dropna() for c in matrix.columns], showfliers=False)
        ax.set_title(title); ax.set_ylabel("Abundance")
    axes[1].set_xticks(range(1, len(before.columns) + 1)); axes[1].set_xticklabels(before.columns, rotation=90)
    fig.tight_layout(); return save_figure(fig, base)


def plot_filtering(summary, base):
    fig, ax = plt.subplots(figsize=(6, 4)); ax.bar(summary["Stage"], summary["Features"], color="#4f46e5")
    for i, value in enumerate(summary["Features"]): ax.text(i, value, f"{int(value):,}", ha="center", va="bottom")
    ax.set_ylabel("Features"); fig.tight_layout(); return save_figure(fig, base)


def plot_qc_rsd(qc_table, base):
    fig, ax = plt.subplots(figsize=(7, 4.5)); data = qc_table["QC RSD"].replace([np.inf, -np.inf], np.nan).dropna()
    ax.hist(data.clip(upper=data.quantile(.99) if len(data) else 1), bins=40, color="#0891b2", edgecolor="white")
    ax.axvline(.2, color="#059669", ls="--", label="20%"); ax.axvline(.3, color="#d97706", ls="--", label="30%")
    ax.set_xlabel("QC RSD"); ax.set_ylabel("Features"); ax.legend(); fig.tight_layout(); return save_figure(fig, base)


def plot_correlation_heatmap(matrix, base, title="Correlation heatmap"):
    fig, ax = plt.subplots(figsize=(max(5, matrix.shape[1] * .45), max(4, matrix.shape[0] * .35)))
    image = ax.imshow(matrix.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(matrix.shape[1])); ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(range(matrix.shape[0])); ax.set_yticklabels(matrix.index); ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Correlation"); fig.tight_layout(); return save_figure(fig, base)


def plot_scores(scores, variance_or_metrics, x_col, y_col, base, qc=True, ellipses=True, labels=False, title=None):
    fig, ax = plt.subplots(figsize=(7, 5.5)); groups = list(pd.unique(scores.get("Group", pd.Series("All", index=scores.index)).fillna("")))
    color_map = {group: GROUP_COLORS[index % len(GROUP_COLORS)] for index, group in enumerate(groups)}
    for group in groups:
        block = scores[scores.get("Group", "All").fillna("").astype(str).eq(str(group))]
        is_qc = str(group).upper() == "QC" or block.get("Sample type", pd.Series(index=block.index, dtype=object)).astype(str).str.lower().eq("qc").all()
        marker = "X" if qc and is_qc else "o"; color = "#111827" if qc and is_qc else color_map[group]
        ax.scatter(block[x_col], block[y_col], label=group or "Unassigned", color=color, marker=marker, s=55, alpha=.9)
        if ellipses and not is_qc: _ellipse(ax, block[x_col].to_numpy(), block[y_col].to_numpy(), color)
        if labels:
            for _, row in block.iterrows(): ax.annotate(str(row["Sample ID"]), (row[x_col], row[y_col]), fontsize=7)
    variance_lookup = {}
    if isinstance(variance_or_metrics, pd.DataFrame) and {"Component", "Explained variance"}.issubset(variance_or_metrics.columns):
        variance_lookup = variance_or_metrics.set_index("Component")["Explained variance"].to_dict()
    x_label = f"{x_col} ({variance_lookup[x_col]:.1%})" if x_col in variance_lookup else x_col
    y_label = f"{y_col} ({variance_lookup[y_col]:.1%})" if y_col in variance_lookup else y_col
    ax.axhline(0, color="#d1d5db", lw=.7); ax.axvline(0, color="#d1d5db", lw=.7)
    ax.set_xlabel(x_label); ax.set_ylabel(y_label); ax.set_title(title or "Score plot"); ax.legend(frameon=False)
    fig.tight_layout(); return save_figure(fig, base)


def plot_pca(pca, output_dir):
    out = Path(output_dir); scores, loadings, variance = pca["scores"], pca["loadings"], pca["variance"]
    paths = {"PCA_PC1_PC2": plot_scores(scores, variance, "PC1", "PC2", out / "Fig_PCA-01_scores_PC1_PC2", title=None)}
    if "PC3" in scores: paths["PCA_PC1_PC3"] = plot_scores(scores, variance, "PC1", "PC3", out / "Fig_PCA-02_scores_PC1_PC3", title=None)
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.bar(variance["Component"], variance["Explained variance"], color="#4f46e5", label="Explained")
    ax.plot(variance["Component"], variance["Cumulative variance"], color="#dc2626", marker="o", label="Cumulative"); ax.set_ylabel("Variance ratio"); ax.legend(); fig.tight_layout()
    paths["PCA_scree"] = save_figure(fig, out / "Fig_PCA-03_scree")
    top = loadings.assign(abs_loading=loadings["PC1"].abs()).nlargest(min(20, len(loadings)), "abs_loading").sort_values("PC1")
    fig, ax = plt.subplots(figsize=(7, max(4, len(top) * .25))); ax.barh(top["Feature ID"], top["PC1"], color=np.where(top["PC1"] >= 0, "#dc2626", "#2563eb")); ax.set_xlabel("PC1 loading"); fig.tight_layout()
    paths["PCA_loadings"] = save_figure(fig, out / "Fig_PCA-04_top_loadings")
    return paths


def plot_plsda(pls, output_dir):
    out = Path(output_dir); paths = {}
    if {"PLS1", "PLS2"}.issubset(pls["scores"].columns): paths["PLS_scores"] = plot_scores(pls["scores"], pls["metrics"], "PLS1", "PLS2", out / "Fig_PLS-01_scores", title=None)
    cv = pls["cv"]; fig, ax1 = plt.subplots(figsize=(7, 4.5)); ax1.plot(cv["Components"], cv["Q2"], marker="o", color="#4f46e5", label="Q2"); ax1.plot(cv["Components"], cv["Balanced accuracy"], marker="s", color="#059669", label="Balanced accuracy"); ax1.axhline(0, color="#d1d5db"); ax1.legend(); ax1.set_xlabel("Components"); fig.tight_layout(); paths["PLS_cv"] = save_figure(fig, out / "Fig_PLS-02_cv")
    perm = pls["permutation"]; true_q2 = float(pls["metrics"].iloc[0]["Q2"]); fig, ax = plt.subplots(figsize=(7, 4.5)); ax.hist(perm["Q2"], bins=30, color="#94a3b8", label="Permuted Q2"); ax.axvline(true_q2, color="#dc2626", lw=2, label=f"Observed Q2={true_q2:.3f}"); ax.legend(); ax.set_xlabel("Q2"); fig.tight_layout(); paths["PLS_permutation"] = save_figure(fig, out / "Fig_PLS-03_permutation")
    top = pls["vip"].head(20).sort_values("VIP"); fig, ax = plt.subplots(figsize=(7, max(4, len(top) * .25))); ax.barh(top["Feature ID"], top["VIP"], color="#7c3aed"); ax.axvline(1, color="#475569", ls="--"); ax.set_xlabel("VIP (interpretive, not significance)"); fig.tight_layout(); paths["PLS_vip"] = save_figure(fig, out / "Fig_PLS-04_VIP")
    return paths


def plot_volcano(table, base, p_column="FDR", exploratory=False, fdr_threshold=.05, fc_threshold=1.2):
    data = table.copy(); y = -np.log10(pd.to_numeric(data[p_column], errors="coerce").clip(lower=1e-300)); x = pd.to_numeric(data["log2FC"], errors="coerce")
    formal = (pd.to_numeric(data[p_column], errors="coerce") < fdr_threshold) & (pd.to_numeric(data["FC"], errors="coerce") >= fc_threshold)
    down = (pd.to_numeric(data[p_column], errors="coerce") < fdr_threshold) & (pd.to_numeric(data["FC"], errors="coerce") <= 1 / fc_threshold)
    colors = np.where(formal, "#dc2626", np.where(down, "#2563eb", "#94a3b8"))
    fig, ax = plt.subplots(figsize=(7, 5.5)); ax.scatter(x, y, c=colors, s=10, alpha=.7, edgecolors="none")
    ax.axvline(np.log2(fc_threshold), color="#475569", ls="--"); ax.axvline(-np.log2(fc_threshold), color="#475569", ls="--"); ax.axhline(-np.log10(fdr_threshold), color="#475569", ls="--")
    ax.set_xlabel("log2FC"); ax.set_ylabel(f"-log10({p_column})"); fig.tight_layout(); return save_figure(fig, base)


def plot_heatmap(abundance, metadata, feature_ids, feature_labels, base, top_n=30):
    ids = [x for x in feature_ids if x in abundance.index][:top_n]; samples = [x for x in metadata[metadata["Sample type"].astype(str).str.lower().ne("qc")]["Sample ID"].astype(str) if x in abundance.columns]
    if not ids or not samples: return {}
    matrix = abundance.loc[ids, samples].apply(pd.to_numeric, errors="coerce"); z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, 1), axis=0).fillna(0)
    if len(z) > 2: z = z.iloc[leaves_list(linkage(z.to_numpy(), method="average", metric="euclidean"))]
    fig, ax = plt.subplots(figsize=(max(7, len(samples) * .35), max(5, len(z) * .22))); image = ax.imshow(z, cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(range(len(samples))); ax.set_xticklabels(samples, rotation=90); labels = [feature_labels.get(i, i) for i in z.index]; ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7); fig.colorbar(image, ax=ax, label="Row Z-score"); fig.tight_layout(); return save_figure(fig, base)


def plot_top_metabolites(abundance, metadata, differential, feature_labels, base, top_n=12):
    data = differential.copy()
    data["rank_q"] = pd.to_numeric(data["FDR"], errors="coerce").fillna(1)
    data["rank_effect"] = pd.to_numeric(data["log2FC"], errors="coerce").abs().fillna(0)
    if "Identified" not in data: data["Identified"] = False
    ids = data.sort_values(["rank_q", "Identified", "rank_effect"], ascending=[True, False, False])["Feature ID"].astype(str).drop_duplicates().head(top_n).tolist()
    bio = metadata[metadata["Sample type"].astype(str).str.lower().ne("qc")]
    groups = [str(value) for value in pd.unique(bio["Group"])]
    columns = 3; rows = int(np.ceil(len(ids) / columns))
    if not ids: return {}
    fig, axes = plt.subplots(rows, columns, figsize=(12, max(3.2, rows * 3)), squeeze=False)
    for ax, feature_id in zip(axes.ravel(), ids):
        blocks = []
        for group in groups:
            sample_ids = bio.loc[bio["Group"].astype(str).eq(group), "Sample ID"].astype(str)
            blocks.append(pd.to_numeric(abundance.loc[feature_id, abundance.columns.intersection(sample_ids)], errors="coerce").dropna().to_numpy())
        ax.boxplot(blocks, tick_labels=groups, showfliers=False)
        for index, values in enumerate(blocks, 1):
            jitter = np.linspace(-.08, .08, len(values)) if len(values) else []
            ax.scatter(index + jitter, values, color=GROUP_COLORS[(index - 1) % len(GROUP_COLORS)], s=14, alpha=.75)
        ax.set_title(feature_labels.get(feature_id, feature_id), fontsize=9)
    for ax in axes.ravel()[len(ids):]: ax.axis("off")
    fig.tight_layout(); return save_figure(fig, base)


def plot_interaction_panels(abundance, metadata, results, factor_a, factor_b, feature_labels, base, top_n=6, fdr_threshold=.05):
    ids = results.loc[pd.to_numeric(results["q_AxB"], errors="coerce") < fdr_threshold].sort_values("q_AxB")["Feature ID"].astype(str).head(top_n).tolist()
    if not ids: return {}
    meta = metadata[metadata["Sample type"].astype(str).str.lower().ne("qc")].set_index("Sample ID")
    a_levels = list(map(str, pd.unique(meta[factor_a]))); b_levels = list(map(str, pd.unique(meta[factor_b])))
    fig, axes = plt.subplots(int(np.ceil(len(ids) / 2)), 2, figsize=(10, max(3.5, np.ceil(len(ids) / 2) * 3.2)), squeeze=False)
    for ax, feature_id in zip(axes.ravel(), ids):
        for index, b in enumerate(b_levels):
            means = []; sems = []
            for a in a_levels:
                sample_ids = meta.index[meta[factor_a].astype(str).eq(a) & meta[factor_b].astype(str).eq(b)].intersection(abundance.columns)
                values = pd.to_numeric(abundance.loc[feature_id, sample_ids], errors="coerce").dropna()
                means.append(values.mean()); sems.append(values.sem())
            ax.errorbar(a_levels, means, yerr=sems, marker="o", capsize=3, label=b, color=GROUP_COLORS[index % len(GROUP_COLORS)])
        ax.set_title(feature_labels.get(feature_id, feature_id), fontsize=9); ax.set_xlabel(factor_a); ax.set_ylabel("Mean ± SEM")
    for ax in axes.ravel()[len(ids):]: ax.axis("off")
    axes.ravel()[0].legend(title=factor_b, frameon=False); fig.tight_layout(); return save_figure(fig, base)


def plot_two_factor_heatmap(results, feature_labels, base, top_n=30):
    data = results.copy()
    q_columns = ["q_A", "q_B", "q_AxB"]
    data["best_q"] = data[q_columns].min(axis=1, skipna=True)
    data = data.sort_values("best_q").head(top_n).set_index("Feature ID")
    if data.empty: return {}
    encoded = pd.DataFrame({label: -np.log10(pd.to_numeric(data[column], errors="coerce").clip(lower=1e-300)) for label, column in zip(["A", "B", "A×B"], q_columns)}, index=data.index)
    fig, ax = plt.subplots(figsize=(6, max(5, len(encoded) * .22))); image = ax.imshow(encoded, cmap="viridis", aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels(encoded.columns); ax.set_yticks(range(len(encoded))); ax.set_yticklabels([feature_labels.get(str(i), str(i)) for i in encoded.index], fontsize=7)
    fig.colorbar(image, ax=ax, label="-log10(FDR)"); fig.tight_layout(); return save_figure(fig, base)


def plot_correlation_scatterplots(abundance, phenotype, correlations, feature_labels, base, top_n=9):
    top = correlations.dropna(subset=["FDR"]).sort_values(["FDR", "P"]).head(top_n)
    if top.empty: return {}
    phenotype = phenotype.copy(); phenotype["Sample ID"] = phenotype["Sample ID"].astype(str); phenotype = phenotype.set_index("Sample ID")
    fig, axes = plt.subplots(int(np.ceil(len(top) / 3)), 3, figsize=(12, max(3.2, np.ceil(len(top) / 3) * 3.1)), squeeze=False)
    for ax, (_, row) in zip(axes.ravel(), top.iterrows()):
        common = phenotype.index.intersection(abundance.columns)
        x = pd.to_numeric(phenotype.loc[common, row["Phenotype"]], errors="coerce"); y = pd.to_numeric(abundance.loc[str(row["Feature ID"]), common], errors="coerce")
        valid = x.notna() & y.notna(); ax.scatter(x[valid], y[valid], color="#4f46e5", s=22)
        if valid.sum() >= 2:
            slope, intercept = np.polyfit(x[valid], y[valid], 1); grid = np.linspace(x[valid].min(), x[valid].max(), 50); ax.plot(grid, slope * grid + intercept, color="#dc2626", lw=1)
        ax.set_xlabel(str(row["Phenotype"])); ax.set_ylabel(feature_labels.get(str(row["Feature ID"]), str(row["Feature ID"])), fontsize=8); ax.set_title(f"r={row['Correlation']:.2f}, q={row['FDR']:.3g}", fontsize=9)
    for ax in axes.ravel()[len(top):]: ax.axis("off")
    fig.tight_layout(); return save_figure(fig, base)


def plot_two_factor_overview(results, effect, base):
    data = results[["Feature ID", f"effect_{effect}", f"q_{effect}"]].dropna().copy(); data = data.sort_values(f"q_{effect}").head(50); data["score"] = -np.log10(data[f"q_{effect}"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(7, 5)); ax.scatter(data[f"effect_{effect}"], data["score"], color="#4f46e5", s=24); ax.set_xlabel("Partial eta-squared"); ax.set_ylabel("-log10(FDR)"); fig.tight_layout(); return save_figure(fig, base)


def plot_annotation(stats_tables, base_dir):
    outputs = {}; base_dir = Path(base_dir)
    for key, table, label in [("coverage", stats_tables["coverage"], "Status"), ("confidence", stats_tables["confidence"], "Annotation level"), ("ion_mode", stats_tables["ion_mode"], "Ion mode")]:
        fig, ax = plt.subplots(figsize=(6, 4)); ax.bar(table[label].astype(str), table["Count"], color=GROUP_COLORS[:len(table)]); ax.tick_params(axis="x", rotation=30); ax.set_ylabel("Features"); fig.tight_layout(); outputs[key] = save_figure(fig, base_dir / f"Fig_ANNO-{['coverage','confidence','ion_mode'].index(key)+1:02d}_{key}")
    return outputs


def plot_pathway(enrichment, base):
    significant = enrichment[enrichment["FDR"] < .05].head(20) if not enrichment.empty else enrichment
    if significant.empty: return {}
    data = significant.sort_values("FDR", ascending=False); fig, ax = plt.subplots(figsize=(8, max(4, len(data) * .3))); scatter = ax.scatter(data["Enrichment ratio"], data["Pathway"], s=30 + data["Hit count"] * 12, c=-np.log10(data["FDR"].clip(lower=1e-300)), cmap="viridis"); fig.colorbar(scatter, ax=ax, label="-log10(FDR)"); ax.set_xlabel("Enrichment ratio"); fig.tight_layout(); return save_figure(fig, base)
