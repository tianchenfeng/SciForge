import os
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

FONT = "Microsoft YaHei"


def parse_args():
    parser = argparse.ArgumentParser(
        description="菌属-环境指标相关分析（总图 + 因素一分图 + 因素二分图）"
    )
    parser.add_argument("--input_dir", required=True, help="输入目录（01_input）")
    parser.add_argument("--output_dir", required=True, help="输出目录（02_output）")
    parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
    return parser.parse_args()


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def find_existing_file(candidates):
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def detect_sample_column(df: pd.DataFrame):
    for col in ["SampleID", "sample-id", "sample_id", "sample", "id"]:
        if col in df.columns:
            return col
    return df.columns[0]


def load_metadata(input_dir: str):
    path = find_existing_file([
        os.path.join(input_dir, "qiime2", "metadata.tsv"),
        os.path.join(input_dir, "metadata.tsv"),
    ])
    if path is None:
        raise FileNotFoundError("找不到 metadata.tsv")
    df = pd.read_csv(path, sep="\t")
    return df, path


def load_env_data(input_dir: str):
    path = find_existing_file([
        os.path.join(input_dir, "data.csv"),
        os.path.join(input_dir, "qiime2", "data.csv"),
    ])
    if path is None:
        raise FileNotFoundError("找不到 data.csv（环境/发酵指标表）")
    df = pd.read_csv(path)
    return df, path


def load_genus_table(output_dir: str):
    path = find_existing_file([
        os.path.join(output_dir, "03_taxa", "Sample_RelAbundance_Genus.tsv"),
        os.path.join(output_dir, "03_taxa", "sample_rel_abundance_genus.tsv"),
    ])
    if path is None:
        raise FileNotFoundError("找不到 Sample_RelAbundance_Genus.tsv，请先运行 03_taxonomy.py")
    df = pd.read_csv(path, sep="\t")
    return df, path


def get_factor_columns(metadata: pd.DataFrame):
    sample_col = detect_sample_column(metadata)
    other_cols = [c for c in metadata.columns if c != sample_col]

    factor1 = other_cols[0] if len(other_cols) >= 1 else None
    factor2 = other_cols[1] if len(other_cols) >= 2 else None
    return sample_col, factor1, factor2


def preprocess_env_table(env_df: pd.DataFrame, sample_col: str):
    env_df = env_df.copy()
    if sample_col not in env_df.columns:
        env_df = env_df.rename(columns={env_df.columns[0]: sample_col})

    env_df[sample_col] = env_df[sample_col].astype(str)

    env_cols = [c for c in env_df.columns if c != sample_col]
    for col in env_cols:
        env_df[col] = pd.to_numeric(env_df[col], errors="coerce")

    env_cols = [c for c in env_cols if pd.api.types.is_numeric_dtype(env_df[c])]
    if len(env_cols) == 0:
        raise ValueError("data.csv 中未找到可用于相关分析的数值型指标列。")

    return env_df, env_cols


def preprocess_genus_table(genus_df: pd.DataFrame, sample_col: str):
    """
    关键修复：
    03_taxa 输出的 Sample_RelAbundance_Genus.tsv 是：
    行 = 菌属
    列 = 样本

    这里统一转换为：
    行 = 样本
    列 = 菌属
    """
    genus_df = genus_df.copy()

    # 第一列是菌属名
    genus_name_col = genus_df.columns[0]

    # 其余列应该是样本列
    sample_like_cols = list(genus_df.columns[1:])
    if len(sample_like_cols) == 0:
        raise ValueError("Sample_RelAbundance_Genus.tsv 不包含样本列。")

    # 转置：行变样本，列变菌属
    genus_wide = genus_df.set_index(genus_name_col).T.reset_index()
    genus_wide = genus_wide.rename(columns={"index": sample_col})

    genus_wide[sample_col] = genus_wide[sample_col].astype(str)

    genus_cols = [c for c in genus_wide.columns if c != sample_col]
    for col in genus_cols:
        genus_wide[col] = pd.to_numeric(genus_wide[col], errors="coerce")

    genus_wide = genus_wide.dropna(axis=0, how="all", subset=genus_cols).fillna(0)

    return genus_wide, genus_cols


def filter_genera(genus_df: pd.DataFrame, genus_cols: list[str], min_mean_abundance=0.01, min_prevalence=0.2):
    keep_cols = []
    n = genus_df.shape[0]
    if n == 0:
        return keep_cols

    for col in genus_cols:
        vals = genus_df[col].values
        mean_abd = np.mean(vals)
        prevalence = np.sum(vals > 0) / n
        if mean_abd >= min_mean_abundance and prevalence >= min_prevalence:
            keep_cols.append(col)

    return keep_cols


def clr_transform(df: pd.DataFrame, pseudo_count=1e-6):
    x = df.astype(float).copy()
    x = x + pseudo_count
    gm = np.exp(np.mean(np.log(x), axis=1))
    clr = np.log(x.div(gm, axis=0))
    return clr


def corr_pvalue_matrix(genus_mat: pd.DataFrame, env_df: pd.DataFrame):
    genus_names = list(genus_mat.columns)
    env_names = list(env_df.columns)

    rho = pd.DataFrame(index=genus_names, columns=env_names, dtype=float)
    pval = pd.DataFrame(index=genus_names, columns=env_names, dtype=float)

    for g in genus_names:
        for e in env_names:
            x = genus_mat[g].values
            y = env_df[e].values

            valid = ~(pd.isna(x) | pd.isna(y))
            if valid.sum() < 3:
                rho.loc[g, e] = np.nan
                pval.loc[g, e] = np.nan
                continue

            r, p = spearmanr(x[valid], y[valid])
            rho.loc[g, e] = r
            pval.loc[g, e] = p

    return rho, pval


def fdr_adjust_matrix(pval_df: pd.DataFrame):
    flat = pval_df.values.flatten()
    mask = ~np.isnan(flat)

    adj = np.full_like(flat, np.nan, dtype=float)
    if mask.sum() > 0:
        _, qvals, _, _ = multipletests(flat[mask], method="fdr_bh")
        adj[mask] = qvals

    qval_df = pd.DataFrame(
        adj.reshape(pval_df.shape),
        index=pval_df.index,
        columns=pval_df.columns
    )
    return qval_df


def significance_star(q):
    if pd.isna(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def build_annotation_df(qval_df: pd.DataFrame):
    ann = qval_df.copy().astype(str)
    for i in qval_df.index:
        for j in qval_df.columns:
            ann.loc[i, j] = significance_star(qval_df.loc[i, j])
    return ann


def sort_genera_by_signal(rho_df: pd.DataFrame):
    score = rho_df.abs().max(axis=1).sort_values(ascending=False)
    return list(score.index)


def plot_heatmap(rho_df: pd.DataFrame, ann_df: pd.DataFrame, title: str, out_png: str, out_tsv: str):
    ordered_rows = sort_genera_by_signal(rho_df)
    rho_df = rho_df.loc[ordered_rows]
    ann_df = ann_df.loc[ordered_rows]

    rho_df.to_csv(out_tsv, sep="\t", encoding="utf-8-sig")

    n_rows, n_cols = rho_df.shape
    fig_w = max(10, n_cols * 1.0 + 2)
    fig_h = max(8, n_rows * 0.42 + 2)

    plt.figure(figsize=(fig_w, fig_h))
    sns.set(style="white")

    ax = sns.heatmap(
        rho_df,
        cmap="RdBu_r",
        center=0,
        vmin=-0.8,
        vmax=0.8,
        linewidths=0.3,
        linecolor="#f0f0f0",
        annot=ann_df,
        fmt="",
        annot_kws={"fontsize": 9, "color": "white", "weight": "bold"},
        cbar_kws={"label": "Spearman Correlation"}
    )

    ax.set_title(title, fontsize=14, pad=12, fontname=FONT)
    ax.set_xlabel("")
    ax.set_ylabel("Genus", fontsize=11, fontname=FONT)

    plt.xticks(rotation=45, ha="right", fontsize=10, fontname=FONT)
    plt.yticks(rotation=0, fontsize=10, fontname=FONT)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def sanitize_tag(x):
    return str(x).replace("/", "_").replace("\\", "_").replace(" ", "_")


def build_master_table(metadata: pd.DataFrame, env_df: pd.DataFrame, genus_df: pd.DataFrame, sample_col: str, factor1: str | None, factor2: str | None):
    metadata_ids = set(metadata[sample_col].astype(str))
    env_ids = set(env_df[sample_col].astype(str))
    genus_ids = set(genus_df[sample_col].astype(str))

    common_ids = metadata_ids & env_ids & genus_ids

    print(f"[INFO] metadata 样本数: {len(metadata_ids)}")
    print(f"[INFO] data.csv 样本数: {len(env_ids)}")
    print(f"[INFO] genus 表样本数: {len(genus_ids)}")
    print(f"[INFO] 三表精确匹配共有样本数: {len(common_ids)}")

    if len(common_ids) < 3:
        raise ValueError("三张表按 SampleID 精确匹配后的共有样本数不足 3，无法进行相关分析。")

    keep_meta_cols = [c for c in [sample_col, factor1, factor2] if c is not None]
    metadata_sub = metadata[metadata[sample_col].isin(common_ids)][keep_meta_cols].copy()
    env_sub = env_df[env_df[sample_col].isin(common_ids)].copy()
    genus_sub = genus_df[genus_df[sample_col].isin(common_ids)].copy()

    merged = pd.merge(metadata_sub, env_sub, on=sample_col, how="inner")
    merged = pd.merge(merged, genus_sub, on=sample_col, how="inner")

    print(f"[INFO] 最终合并样本数: {merged.shape[0]}")
    return merged


def analyze_one_subset(sub_df: pd.DataFrame, sample_col: str, env_cols: list[str], factor_cols: list[str], out_png: str, out_tsv: str, title: str):
    exclude_cols = [sample_col] + env_cols + factor_cols
    genus_cols = [c for c in sub_df.columns if c not in exclude_cols]

    if sub_df.shape[0] < 3:
        return False, "样本数不足（<3）"

    keep_cols = filter_genera(
        pd.concat([sub_df[[sample_col]], sub_df[genus_cols]], axis=1),
        genus_cols,
        min_mean_abundance=0.01,
        min_prevalence=0.2
    )
    if len(keep_cols) == 0:
        return False, "筛选后无可用菌属"

    genus_sub = sub_df[keep_cols].copy()
    env_sub = sub_df[env_cols].copy()

    genus_clr = clr_transform(genus_sub)

    rho_df, p_df = corr_pvalue_matrix(genus_clr, env_sub)
    q_df = fdr_adjust_matrix(p_df)
    ann_df = build_annotation_df(q_df)

    plot_heatmap(rho_df, ann_df, title, out_png, out_tsv)
    return True, f"n_samples={sub_df.shape[0]}, n_genera={len(keep_cols)}"


def main():
    args = parse_args()

    out_dir = os.path.join(args.output_dir, "07_env_corr")
    ensure_dir(out_dir)

    metadata, metadata_path = load_metadata(args.input_dir)
    env_df, env_path = load_env_data(args.input_dir)
    genus_df_raw, genus_path = load_genus_table(args.output_dir)

    sample_col, factor1, factor2 = get_factor_columns(metadata)

    metadata[sample_col] = metadata[sample_col].astype(str)
    env_df, env_cols = preprocess_env_table(env_df, sample_col)
    genus_df, genus_cols = preprocess_genus_table(genus_df_raw, sample_col)

    merged = build_master_table(metadata, env_df, genus_df, sample_col, factor1, factor2)

    factor_cols = [c for c in [factor1, factor2] if c is not None]
    summary_records = []

    # 总图
    overall_png = os.path.join(out_dir, "env_corr_overall.png")
    overall_tsv = os.path.join(out_dir, "env_corr_overall.tsv")
    ok, msg = analyze_one_subset(
        sub_df=merged,
        sample_col=sample_col,
        env_cols=env_cols,
        factor_cols=factor_cols,
        out_png=overall_png,
        out_tsv=overall_tsv,
        title="Genus vs Environment Correlation (Overall)"
    )
    summary_records.append({
        "group_type": "overall",
        "group_name": "overall",
        "status": "OK" if ok else "SKIP",
        "detail": msg,
        "output_png": os.path.basename(overall_png),
        "output_tsv": os.path.basename(overall_tsv),
    })

    # 因素一分图
    if factor1 is not None:
        factor1_dir = os.path.join(out_dir, f"by_{factor1}")
        ensure_dir(factor1_dir)

        for lv in merged[factor1].dropna().unique():
            sub = merged[merged[factor1] == lv].copy()
            tag = sanitize_tag(lv)

            out_png = os.path.join(factor1_dir, f"env_corr_{factor1}_{tag}.png")
            out_tsv = os.path.join(factor1_dir, f"env_corr_{factor1}_{tag}.tsv")

            ok, msg = analyze_one_subset(
                sub_df=sub,
                sample_col=sample_col,
                env_cols=env_cols,
                factor_cols=factor_cols,
                out_png=out_png,
                out_tsv=out_tsv,
                title=f"Genus vs Environment Correlation ({factor1}={lv})"
            )

            summary_records.append({
                "group_type": factor1,
                "group_name": str(lv),
                "status": "OK" if ok else "SKIP",
                "detail": msg,
                "output_png": os.path.relpath(out_png, out_dir),
                "output_tsv": os.path.relpath(out_tsv, out_dir),
            })

    # 因素二分图
    if factor2 is not None:
        factor2_dir = os.path.join(out_dir, f"by_{factor2}")
        ensure_dir(factor2_dir)

        for lv in merged[factor2].dropna().unique():
            sub = merged[merged[factor2] == lv].copy()
            tag = sanitize_tag(lv)

            out_png = os.path.join(factor2_dir, f"env_corr_{factor2}_{tag}.png")
            out_tsv = os.path.join(factor2_dir, f"env_corr_{factor2}_{tag}.tsv")

            ok, msg = analyze_one_subset(
                sub_df=sub,
                sample_col=sample_col,
                env_cols=env_cols,
                factor_cols=factor_cols,
                out_png=out_png,
                out_tsv=out_tsv,
                title=f"Genus vs Environment Correlation ({factor2}={lv})"
            )

            summary_records.append({
                "group_type": factor2,
                "group_name": str(lv),
                "status": "OK" if ok else "SKIP",
                "detail": msg,
                "output_png": os.path.relpath(out_png, out_dir),
                "output_tsv": os.path.relpath(out_tsv, out_dir),
            })

    pd.DataFrame(summary_records).to_csv(
        os.path.join(out_dir, "env_corr_summary.tsv"),
        sep="\t",
        index=False,
        encoding="utf-8-sig"
    )

    print("[OK] 菌属-环境指标相关分析完成")
    print(f"  -> metadata: {metadata_path}")
    print(f"  -> env data: {env_path}")
    print(f"  -> genus abundance: {genus_path}")
    print(f"  -> factor1: {factor1}")
    print(f"  -> factor2: {factor2}")
    print(f"  -> env columns: {', '.join(env_cols)}")
    print(f"  -> output: {out_dir}")


if __name__ == "__main__":
    main()