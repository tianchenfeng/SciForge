import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="06_genus_function_association")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_ROOT = args.output_dir
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "06_assoc")
ANALYSIS_MODE = args.mode

os.makedirs(OUTPUT_DIR, exist_ok=True)

META_PATH_CANDIDATES = [
    os.path.join(INPUT_DIR, "qiime2", "metadata.tsv"),
    os.path.join(INPUT_DIR, "qiime2", "metadata.csv"),
    os.path.join(INPUT_DIR, "metadata.tsv"),
    os.path.join(INPUT_DIR, "metadata.csv"),
]

# 直接读取前面步骤在“当前输出目录”下生成的结果
GENUS_TABLE_CANDIDATES = [
    os.path.join(OUTPUT_ROOT, "03_taxa", "Sample_RelAbundance_Genus.tsv"),
]

KO_TABLE_CANDIDATES = [
    os.path.join(OUTPUT_ROOT, "05_picrust2", "KO_InferredMetagenome_CPM.tsv"),
]
EC_TABLE_CANDIDATES = [
    os.path.join(OUTPUT_ROOT, "05_picrust2", "EC_InferredMetagenome_CPM.tsv"),
]

KO_HEATMAP_TOP_CANDIDATES = [
    os.path.join(OUTPUT_ROOT, "05_picrust2", "KO_TopVariableFunctions_GroupMean.tsv"),
]
EC_HEATMAP_TOP_CANDIDATES = [
    os.path.join(OUTPUT_ROOT, "05_picrust2", "EC_TopVariableFunctions_GroupMean.tsv"),
]

FONT = "Arial"

TOP_N_GENUS = 20
TOP_N_FUNCTION = 30

CORR_THRESHOLD = 0.60
FDR_THRESHOLD = 0.05

HEATMAP_FIG_W = 12
HEATMAP_MIN_H = 6
HEATMAP_MAX_H = 16


# ================= 2. 工具函数 =================
def safe_print(msg):
    try:
        print(str(msg), flush=True)
    except UnicodeEncodeError:
        safe_msg = str(msg).encode("gbk", errors="replace").decode("gbk", errors="replace")
        print(safe_msg, flush=True)

def find_first_existing(paths, allow_none=False):
    for p in paths:
        if os.path.exists(p):
            return p
    if allow_none:
        return None
    raise FileNotFoundError("以下候选路径均不存在：\n" + "\n".join(paths))


def detect_sep(file_path):
    with open(file_path, "r", encoding="utf-8-sig") as f:
        first_line = f.readline()
    return "\t" if "\t" in first_line else ","


def clean_text(x):
    if pd.isna(x):
        return ""
    x = str(x)
    x = x.replace("\ufeff", "")
    x = x.replace("\u200b", "")
    x = x.replace("\xa0", " ")
    return x.strip()


def clean_header_text(x):
    x = clean_text(x)
    x = re.sub(r"^[\[\]\(\)]+", "", x)
    x = re.sub(r"[\[\]\(\)]+$", "", x)
    return x.strip()


def read_table_auto(path, index_col=0):
    sep = detect_sep(path)
    df = pd.read_csv(path, sep=sep, index_col=index_col, encoding="utf-8-sig")
    df.index = df.index.astype(str).map(clean_text)
    df.columns = [clean_header_text(c) for c in df.columns]
    return df


def natural_sort(values):
    vals = [clean_text(v) for v in pd.unique(pd.Series(values).astype(str))]

    def try_float(x):
        try:
            return float(x)
        except Exception:
            return None

    parsed = [try_float(v) for v in vals]
    if all(v is not None for v in parsed):
        return [str(x) for x in sorted(vals, key=lambda z: float(z))]

    def natural_key(s):
        return [int(text) if text.isdigit() else text.lower()
                for text in re.split(r"(\d+)", str(s))]

    return sorted(vals, key=natural_key)


def order_row_factor(values):
    vals = [clean_text(v) for v in pd.unique(pd.Series(values).astype(str))]
    priority = []
    for x in ["C", "CON", "Control", "N"]:
        if x in vals:
            priority.append(x)
    rest = [v for v in vals if v not in priority]
    return priority + natural_sort(rest)


def build_sample_order_from_meta(meta, factor1, factor2):
    return meta.index.tolist()


def shorten_name(name, max_len=55):
    name = clean_text(name)
    if len(name) <= max_len:
        return name
    return name[:max_len] + "..."


# ================= 3. metadata =================
def read_metadata_auto(meta_path):
    sep = detect_sep(meta_path)
    meta = pd.read_csv(meta_path, sep=sep, dtype=str, encoding="utf-8-sig")

    if meta.shape[0] > 0 and str(meta.iloc[0, 0]).startswith("#q2"):
        meta = meta.iloc[1:].copy()

    meta.columns = [clean_header_text(c) for c in meta.columns]

    min_columns = 2 if ANALYSIS_MODE == "single_factor" else 3
    if meta.shape[1] < min_columns:
        expected = "样本名 + 一个因素" if ANALYSIS_MODE == "single_factor" else "样本名 + 两个因素"
        raise ValueError(
            f"metadata 至少需要 {min_columns} 列：{expected}。\n当前列名：{list(meta.columns)}"
        )

    sample_col = meta.columns[0]
    factor1 = clean_header_text(meta.columns[1])
    factor2 = None if ANALYSIS_MODE == "single_factor" else clean_header_text(meta.columns[2])

    factor_columns = [factor1] + ([factor2] if factor2 else [])
    meta = meta[[sample_col] + factor_columns].copy()
    meta = meta.rename(columns={sample_col: "SampleID"})

    meta["SampleID"] = meta["SampleID"].astype(str).map(clean_text)
    meta[factor1] = meta[factor1].astype(str).map(clean_text)
    if factor2:
        meta[factor2] = meta[factor2].astype(str).map(clean_text)

    meta = meta.loc[~meta["SampleID"].duplicated(keep="first")].copy()
    meta = meta.set_index("SampleID")

    return meta, factor1, factor2


# ================= 4. 数据准备 =================
def load_base_tables():
    safe_print("[Step 1] Loading metadata / genus tables ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    genus_path = find_first_existing(GENUS_TABLE_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> genus table: {genus_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    genus_df = read_table_auto(genus_path, index_col=0)

    valid_samples = [s for s in build_sample_order_from_meta(meta, factor1, factor2) if s in genus_df.columns]
    if len(valid_samples) == 0:
        raise ValueError("metadata 与 genus 表没有共同样本。")

    meta = meta.loc[valid_samples].copy()
    genus_df = genus_df.loc[:, valid_samples].copy()

    safe_print(f"  -> Common samples: {len(valid_samples)}")
    safe_print(f"  -> Row factor (第2列): {factor1}")
    if factor2:
        safe_print(f"  -> Col factor (第3列): {factor2}")
    else:
        safe_print("  -> Analysis mode: single factor")

    return meta, genus_df, factor1, factor2


def load_function_table(prefix):
    if prefix == "KO":
        func_path = find_first_existing(KO_TABLE_CANDIDATES, allow_none=True)
        top_path = find_first_existing(KO_HEATMAP_TOP_CANDIDATES, allow_none=True)
    elif prefix == "EC":
        func_path = find_first_existing(EC_TABLE_CANDIDATES, allow_none=True)
        top_path = find_first_existing(EC_HEATMAP_TOP_CANDIDATES, allow_none=True)
    else:
        raise ValueError(f"Unsupported prefix: {prefix}")

    if func_path is None:
        safe_print(f"  [Skip] {prefix} 功能表不存在。")
        return None, None

    safe_print(f"  -> {prefix} function table: {func_path}")
    safe_print(f"  -> {prefix} top-function reference: {top_path if top_path else 'Not found'}")

    func_df = read_table_auto(func_path, index_col=0)
    return func_df, top_path


def select_top_genera(genus_df, top_n=TOP_N_GENUS):
    mean_abd = genus_df.mean(axis=1).sort_values(ascending=False)
    top_genera = mean_abd.head(min(top_n, len(mean_abd))).index.tolist()
    return genus_df.loc[top_genera].copy()


def select_top_functions(func_df, top_path=None, top_n=TOP_N_FUNCTION):
    """
    func_df: sample × function
    优先依据 05 的 top variable function groupmean 表选功能
    """
    if top_path is not None and os.path.exists(top_path):
        top_group_df = read_table_auto(top_path, index_col=0)
        # top_group_df 是 function × groups
        top_funcs = top_group_df.index.astype(str).map(clean_text).tolist()
        matched = [f for f in top_funcs if f in func_df.columns]
        if len(matched) > 0:
            return func_df.loc[:, matched].copy()

    func_sd = func_df.std(axis=0).sort_values(ascending=False)
    top_funcs = func_sd.head(min(top_n, len(func_sd))).index.tolist()
    return func_df.loc[:, top_funcs].copy()


# ================= 5. 相关分析 =================
def compute_spearman_matrix(genus_df, func_df):
    """
    genus_df: genus × sample
    func_df : sample × function
    """
    safe_print("[Step 2] Computing Spearman correlations ...")

    genera = genus_df.index.tolist()
    functions = func_df.columns.tolist()

    rho_mat = pd.DataFrame(index=genera, columns=functions, dtype=float)
    p_mat = pd.DataFrame(index=genera, columns=functions, dtype=float)

    records = []

    for g in genera:
        x = pd.to_numeric(genus_df.loc[g], errors="coerce").values

        for f in functions:
            y = pd.to_numeric(func_df[f], errors="coerce").values

            valid = np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 3:
                rho, p = np.nan, np.nan
            else:
                try:
                    rho, p = spearmanr(x[valid], y[valid])
                except Exception:
                    rho, p = np.nan, np.nan

            rho_mat.loc[g, f] = rho
            p_mat.loc[g, f] = p

            records.append({
                "Genus": g,
                "Function": f,
                "Spearman_rho": rho,
                "P_value": p
            })

    long_df = pd.DataFrame(records)

    valid_mask = long_df["P_value"].notna()
    fdr_vals = np.full(len(long_df), np.nan)

    if valid_mask.sum() > 0:
        _, fdr_corrected, _, _ = multipletests(
            long_df.loc[valid_mask, "P_value"].values,
            alpha=FDR_THRESHOLD,
            method="fdr_bh"
        )
        fdr_vals[valid_mask.values] = fdr_corrected

    long_df["FDR"] = fdr_vals

    fdr_mat = pd.DataFrame(index=genera, columns=functions, dtype=float)
    for _, row in long_df.iterrows():
        fdr_mat.loc[row["Genus"], row["Function"]] = row["FDR"]

    return rho_mat, p_mat, fdr_mat, long_df


def filter_assoc_results(rho_mat, fdr_mat, corr_thr=CORR_THRESHOLD, fdr_thr=FDR_THRESHOLD):
    mask = (rho_mat.abs() >= corr_thr) & (fdr_mat < fdr_thr)
    filt_rho = rho_mat.where(mask)
    filt_rho = filt_rho.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return filt_rho


def build_star_matrix(rho_mat, fdr_mat):
    star_mat = pd.DataFrame("", index=rho_mat.index, columns=rho_mat.columns, dtype=object)

    for i in rho_mat.index:
        for j in rho_mat.columns:
            p = fdr_mat.loc[i, j]
            if pd.isna(p):
                star = ""
            elif p < 0.001:
                star = "***"
            elif p < 0.01:
                star = "**"
            elif p < 0.05:
                star = "*"
            else:
                star = ""
            star_mat.loc[i, j] = star

    return star_mat


# ================= 6. 导出表 =================
def export_assoc_tables(prefix, rho_mat, p_mat, fdr_mat, long_df, filt_rho):
    safe_print(f"[Step 3] Exporting association tables for {prefix} ...")

    rho_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_SpearmanRho.tsv")
    p_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_Pvalue.tsv")
    fdr_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_FDR.tsv")
    long_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_Association_Long.tsv")
    filt_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_FilteredRho.tsv")

    rho_mat.to_csv(rho_path, sep="\t")
    p_mat.to_csv(p_path, sep="\t")
    fdr_mat.to_csv(fdr_path, sep="\t")
    long_df.to_csv(long_path, sep="\t", index=False)
    filt_rho.to_csv(filt_path, sep="\t")

    safe_print(f"  [OK] Saved: {rho_path}")
    safe_print(f"  [OK] Saved: {p_path}")
    safe_print(f"  [OK] Saved: {fdr_path}")
    safe_print(f"  [OK] Saved: {long_path}")
    safe_print(f"  [OK] Saved: {filt_path}")


# ================= 7. 作图 =================
def plot_assoc_heatmap(prefix, filt_rho, fdr_mat):
    safe_print(f"[Step 4] Plotting {prefix} association heatmap ...")

    if filt_rho.empty:
        safe_print(f"  [Warning] {prefix} 经过 |rho| 和 FDR 筛选后没有保留项，跳过热图。")
        return

    plot_df = filt_rho.copy()
    plot_df.index = [shorten_name(x, 35) for x in plot_df.index]
    plot_df.columns = [shorten_name(x, 45) for x in plot_df.columns]

    star_sub = build_star_matrix(
        rho_mat=filt_rho,
        fdr_mat=fdr_mat.loc[filt_rho.index, filt_rho.columns]
    )
    star_sub.index = plot_df.index
    star_sub.columns = plot_df.columns

    fig_h = min(max(HEATMAP_MIN_H, plot_df.shape[0] * 0.45), HEATMAP_MAX_H)
    fig_w = max(HEATMAP_FIG_W, plot_df.shape[1] * 0.45)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    hm = sns.heatmap(
        plot_df,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        annot=star_sub,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={
            "label": "Spearman rho",
            "shrink": 0.90,
            "pad": 0.03
        },
        ax=ax
    )

    ax.set_title(
        f"{prefix} genus-function associations (|rho| ≥ {CORR_THRESHOLD}, FDR < {FDR_THRESHOLD})",
        fontname=FONT,
        fontsize=13
    )
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.xticks(rotation=45, ha="right", fontsize=9, fontname=FONT)
    plt.yticks(rotation=0, fontsize=9, fontname=FONT)

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10, pad=2)

    fig.subplots_adjust(left=0.28, right=0.90, bottom=0.18, top=0.92)

    png_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_Association_Heatmap.png")
    pdf_path = os.path.join(OUTPUT_DIR, f"{prefix}_Genus_Function_Association_Heatmap.pdf")

    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {png_path}")


# ================= 8. 单套流程 =================
def run_assoc_for_prefix(prefix, meta, genus_df):
    safe_print(f"\n========== Processing {prefix} ==========")

    func_df, top_path = load_function_table(prefix)
    if func_df is None:
        return

    valid_samples = [s for s in meta.index if (s in genus_df.columns) and (s in func_df.index)]
    if len(valid_samples) == 0:
        safe_print(f"  [Skip] {prefix} 与 genus 表无共同样本。")
        return

    meta_sub = meta.loc[valid_samples].copy()
    genus_sub = genus_df.loc[:, valid_samples].copy()
    func_sub = func_df.loc[valid_samples, :].copy()

    genus_top = select_top_genera(genus_sub, top_n=TOP_N_GENUS)
    func_top = select_top_functions(func_sub, top_path=top_path, top_n=TOP_N_FUNCTION)

    genus_top_path = os.path.join(OUTPUT_DIR, f"{prefix}_Top_Genus_Input.tsv")
    func_top_path = os.path.join(OUTPUT_DIR, f"{prefix}_Top_Function_Input.tsv")
    genus_top.to_csv(genus_top_path, sep="\t")
    func_top.to_csv(func_top_path, sep="\t")

    safe_print(f"  [OK] Saved: {genus_top_path}")
    safe_print(f"  [OK] Saved: {func_top_path}")

    rho_mat, p_mat, fdr_mat, long_df = compute_spearman_matrix(genus_top, func_top)
    filt_rho = filter_assoc_results(
        rho_mat,
        fdr_mat,
        corr_thr=CORR_THRESHOLD,
        fdr_thr=FDR_THRESHOLD
    )

    export_assoc_tables(prefix, rho_mat, p_mat, fdr_mat, long_df, filt_rho)
    plot_assoc_heatmap(prefix, filt_rho, fdr_mat)


# ================= 9. 主程序 =================
if __name__ == "__main__":
    meta, genus_df, factor1, factor2 = load_base_tables()

    run_assoc_for_prefix("KO", meta, genus_df)
    run_assoc_for_prefix("EC", meta, genus_df)

    safe_print("\n[OK] [06_assoc] Done!")
