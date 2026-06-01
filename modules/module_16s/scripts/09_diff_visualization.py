import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests


# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="09_stratified_differential_visualization")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "09_diff_viz")
ANALYSIS_MODE = args.mode

os.makedirs(OUTPUT_DIR, exist_ok=True)

META_PATH_CANDIDATES = [
    os.path.join(INPUT_DIR, "qiime2", "metadata.tsv"),
    os.path.join(INPUT_DIR, "qiime2", "metadata.csv"),
    os.path.join(INPUT_DIR, "metadata.tsv"),
    os.path.join(INPUT_DIR, "metadata.csv"),
]

TABLE_PATH_CANDIDATES = [
    os.path.join(INPUT_DIR, "qiime2", "feature-table.tsv"),
    os.path.join(INPUT_DIR, "feature-table.tsv"),
]

TAX_PATH_CANDIDATES = [
    os.path.join(INPUT_DIR, "qiime2", "taxonomy.tsv"),
    os.path.join(INPUT_DIR, "taxonomy.tsv"),
]

FONT = "Arial"
NPG_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"
]

P_THRESHOLD = 0.05
FDR_THRESHOLD = 0.05
LDA_THRESHOLD = 1.5
TOP_N_PLOT = 15
MIN_MEAN_ABUNDANCE = 0.01
MIN_GROUP_SIZE = 2


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

EXCLUDE_TAXA_KEYWORDS = (
    "unclassified",
    "mitochondria",
    "mitochondrion",
    "chloroplast",
    "uncultured",
    "unidentified",
    "unknown",
    "metagenome",
    "ambiguous_taxa",
)


def is_excluded_taxon_name(name):
    s = clean_text(name).lower()
    return any(k in s for k in EXCLUDE_TAXA_KEYWORDS)


def filter_excluded_taxa(df):
    keep_index = [not is_excluded_taxon_name(x) for x in df.index]
    return df.loc[keep_index].copy()

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


def order_factor1(values):
    vals = [clean_text(v) for v in pd.unique(pd.Series(values).astype(str))]
    priority = []
    for x in ["C", "CON", "Control", "N"]:
        if x in vals:
            priority.append(x)
    rest = [v for v in vals if v not in priority]
    return priority + natural_sort(rest)


def shorten_name(name, max_len=45):
    name = clean_text(name)
    if len(name) <= max_len:
        return name
    return name[:max_len] + "..."


def compute_lda_like_score(group_means):
    diff = float(group_means.max() - group_means.min())
    score = np.log10(1 + diff) * 4
    return score


# ================= 3. 数据读取 =================
def read_metadata_auto(meta_path):
    sep = detect_sep(meta_path)
    meta = pd.read_csv(meta_path, sep=sep, dtype=str, encoding="utf-8-sig")

    if meta.shape[0] > 0 and str(meta.iloc[0, 0]).startswith("#q2"):
        meta = meta.iloc[1:].copy()

    meta.columns = [clean_header_text(c) for c in meta.columns]

    if meta.shape[1] < 3:
        raise ValueError(
            f"metadata 至少需要 3 列：样本名 + 两个因素。实际列名：{list(meta.columns)}"
        )

    sample_col = meta.columns[0]
    factor1 = clean_header_text(meta.columns[1])
    factor2 = clean_header_text(meta.columns[2])

    meta = meta[[sample_col, factor1, factor2]].copy()
    meta = meta.rename(columns={sample_col: "SampleID"})

    meta["SampleID"] = meta["SampleID"].astype(str).map(clean_text)
    meta[factor1] = meta[factor1].astype(str).map(clean_text)
    meta[factor2] = meta[factor2].astype(str).map(clean_text)

    meta = meta.loc[~meta["SampleID"].duplicated(keep="first")].copy()
    meta = meta.set_index("SampleID")

    return meta, factor1, factor2


def read_feature_table(table_path):
    ft = pd.read_csv(table_path, sep="\t", skiprows=1, index_col=0)
    ft.index = ft.index.astype(str).map(clean_text)
    ft.columns = [clean_text(c) for c in ft.columns]
    ft = ft.apply(pd.to_numeric, errors="coerce").fillna(0)
    ft = ft.loc[~ft.index.astype(str).str.startswith("#"), :]
    return ft


def read_taxonomy(tax_path):
    tax = pd.read_csv(tax_path, sep="\t", index_col=0)
    tax.index = tax.index.astype(str).map(clean_text)
    tax.columns = [clean_header_text(c) for c in tax.columns]

    if "Taxon" not in tax.columns:
        raise ValueError("taxonomy.tsv 缺少 Taxon 列。")

    tax["Taxon"] = tax["Taxon"].astype(str).map(clean_text)
    tax = tax.loc[~tax.index.duplicated(keep="first")].copy()
    return tax


def parse_genus(taxon_str):
    if pd.isna(taxon_str):
        return "Unclassified"

    parts = [clean_text(x) for x in str(taxon_str).split(";")]
    genus = "Unclassified"

    for part in parts:
        if part.startswith("g__"):
            genus = part.replace("g__", "").strip()
            break

    if genus in ["", "Unclassified"]:
        for item in reversed(parts):
            tmp = item
            for prefix in ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]:
                tmp = tmp.replace(prefix, "")
            tmp = tmp.strip()
            if tmp != "":
                genus = tmp
                break

    bad = {"", "uncultured", "unidentified", "unknown", "metagenome", "ambiguous_taxa"}
    if genus.lower() in bad:
        return "Unclassified"

    return clean_text(genus)


def load_and_clean_data():
    safe_print("[Step 1] Loading & cleaning data ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)
    tax_path = find_first_existing(TAX_PATH_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")
    safe_print(f"  -> taxonomy: {tax_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)
    tax = read_taxonomy(tax_path)

    valid_samples = [s for s in meta.index if s in ft.columns]
    if len(valid_samples) == 0:
        raise ValueError("metadata 与 feature-table 没有共同样本。")

    meta = meta.loc[valid_samples].copy()
    ft = ft.loc[:, valid_samples].copy()

    common_features = ft.index.intersection(tax.index)
    if len(common_features) == 0:
        raise ValueError("feature-table 与 taxonomy 没有共同 Feature ID。")

    ft = ft.loc[common_features].copy()
    tax = tax.loc[common_features].copy()

    # 先从 feature 层面过滤非目标序列
    keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)

    ft_use = ft.loc[keep_features].copy()
    tax_use = tax.loc[keep_features].copy()

    if ft_use.empty:
        raise ValueError("过滤线粒体/叶绿体/未分类后没有有效 feature。")

    tax_use["Genus_Clean"] = tax_use["Taxon"].apply(parse_genus)

    # 先聚合到 genus
    ft_genus = ft_use.groupby(tax_use["Genus_Clean"]).sum()

    # 再保险：聚合后的 genus 名继续关键词过滤
    ft_genus = filter_excluded_taxa(ft_genus)

    if ft_genus.empty:
        raise ValueError("过滤后没有有效 genus。")

    # 关键：过滤后重新计算相对丰度，使每个样本总和为 100%
    ft_genus = ft_genus.div(ft_genus.sum(axis=0).replace(0, np.nan), axis=1) * 100
    ft_genus = ft_genus.fillna(0)

    ft_genus = ft_genus.loc[ft_genus.mean(axis=1) > MIN_MEAN_ABUNDANCE].copy()

    safe_print(f"  -> Samples: {len(valid_samples)}")
    safe_print(f"  -> Retained genera: {ft_genus.shape[0]}")
    safe_print(f"  -> Factor1: {factor1}")
    safe_print(f"  -> Factor2: {factor2}")

    return ft_genus, meta, factor1, factor2


# ================= 4. 统计核心 =================
def run_stat_analysis(df_rel, meta_sub, compare_col):
    groups = natural_sort(meta_sub[compare_col].unique())
    records = []

    if len(groups) < 2:
        return pd.DataFrame()

    raw_pvals = []
    tmp_rows = []

    for taxon in df_rel.index:
        values_by_group = []
        group_means = {}
        valid_group_count = 0

        for g in groups:
            vals = df_rel.loc[taxon, meta_sub.index[meta_sub[compare_col] == g]].values.astype(float)
            vals = vals[np.isfinite(vals)]

            if len(vals) >= MIN_GROUP_SIZE:
                values_by_group.append(vals)
                valid_group_count += 1

            group_means[g] = np.mean(vals) if len(vals) > 0 else np.nan

        if valid_group_count < 2:
            continue

        if np.nansum(list(group_means.values())) <= 0:
            continue

        try:
            _, p = kruskal(*values_by_group)
        except Exception:
            p = np.nan

        means_s = pd.Series(group_means).dropna()
        if means_s.empty:
            continue

        enriched_group = means_s.idxmax()
        lda_like = compute_lda_like_score(means_s)

        row = {
            "Taxon": taxon,
            "P_value": p,
            "Enriched": enriched_group,
            "LDA_like": lda_like,
            **{f"Mean_{g}": group_means[g] for g in groups}
        }

        tmp_rows.append(row)
        raw_pvals.append(p)

    if len(tmp_rows) == 0:
        return pd.DataFrame()

    df = pd.DataFrame(tmp_rows)

    valid_mask = df["P_value"].notna()
    fdr_vals = np.full(len(df), np.nan)
    if valid_mask.sum() > 0:
        _, fdr_corrected, _, _ = multipletests(
            df.loc[valid_mask, "P_value"].values,
            alpha=FDR_THRESHOLD,
            method="fdr_bh"
        )
        fdr_vals[valid_mask.values] = fdr_corrected

    df["FDR"] = fdr_vals
    df = df.sort_values(["FDR", "LDA_like"], ascending=[True, False]).reset_index(drop=True)

    return df


def filter_stat_results(df):
    if df.empty:
        return df
    return df[
        (df["P_value"] < P_THRESHOLD) &
        (df["FDR"] < FDR_THRESHOLD) &
        (df["LDA_like"] >= LDA_THRESHOLD)
    ].copy()


# ================= 5. 绘图 =================
def plot_barplot(df_filt, df_all, compare_col, title_prefix, out_dir):
    if df_filt.empty:
        safe_print(f"  [Info] {title_prefix}: no strict-significant candidates, using fallback top taxa.")
        plot_df = df_all.sort_values(["LDA_like", "P_value"], ascending=[False, True]).head(TOP_N_PLOT).copy()
        suffix = "_Fallback"
    else:
        plot_df = df_filt.sort_values("LDA_like", ascending=True).tail(TOP_N_PLOT).copy()
        suffix = ""

    if plot_df.empty:
        return

    groups = natural_sort(plot_df["Enriched"].unique())
    color_map = dict(zip(groups, (NPG_COLORS * 20)[:len(groups)]))
    colors = [color_map[g] for g in plot_df["Enriched"]]

    fig_h = max(6, len(plot_df) * 0.35)
    fig, ax = plt.subplots(figsize=(9, fig_h))

    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["LDA_like"], color=colors, edgecolor="black", linewidth=0.4)

    ax.set_yticks(y)
    ax.set_yticklabels([shorten_name(x, 45) for x in plot_df["Taxon"]], fontname=FONT, fontsize=9)
    ax.set_xlabel("LDA-like score", fontname=FONT, fontsize=12)
    ax.set_ylabel("")
    ax.set_title(f"{title_prefix}\n(Kruskal + FDR + LDA-like)", fontname=FONT, fontsize=13)

    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[g]) for g in groups]
    ax.legend(handles, groups, title=f"Enriched in {compare_col}", frameon=False,
              bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()

    safe_title = title_prefix.replace(" ", "_").replace(":", "").replace("/", "_")
    png_path = os.path.join(out_dir, f"{safe_title}_Bar{suffix}.png")
    pdf_path = os.path.join(out_dir, f"{safe_title}_Bar{suffix}.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {png_path}")


def plot_heatmap(df_filt, df_all, ft_sub, meta_sub, compare_col, title_prefix, out_dir):
    if df_filt.empty:
        safe_print(f"  [Info] {title_prefix}: no strict-significant candidates, using fallback top taxa heatmap.")
        sig_taxa = df_all.sort_values(["LDA_like", "P_value"], ascending=[False, True]).head(TOP_N_PLOT)["Taxon"].tolist()
        suffix = "_Fallback"
    else:
        sig_taxa = df_filt.sort_values("LDA_like", ascending=False).head(TOP_N_PLOT)["Taxon"].tolist()
        suffix = ""

    if len(sig_taxa) == 0:
        return

    heatmap_data = ft_sub.loc[sig_taxa].T.copy()
    heatmap_data["Group"] = meta_sub[compare_col]

    heatmap_mean = heatmap_data.groupby("Group").mean().T
    group_order = natural_sort(meta_sub[compare_col].unique())
    heatmap_mean = heatmap_mean.loc[:, [g for g in group_order if g in heatmap_mean.columns]]

    def row_zscore(x):
        sd = x.std(ddof=1)
        if pd.isna(sd) or np.isclose(sd, 0, atol=1e-15):
            return pd.Series(np.zeros(len(x)), index=x.index)
        return (x - x.mean()) / sd

    heatmap_z = heatmap_mean.apply(row_zscore, axis=1).clip(-2, 2)

    fig_h = max(6, len(sig_taxa) * 0.35)
    fig, ax = plt.subplots(figsize=(6.5, fig_h))

    hm = sns.heatmap(
        heatmap_z,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        annot=True,
        fmt=".1f",
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "Z-score", "shrink": 0.90, "pad": 0.03},
        ax=ax
    )

    ax.set_title(f"{title_prefix} heatmap", fontname=FONT, fontsize=13)
    ax.set_xlabel(compare_col, fontname=FONT, fontsize=11)
    ax.set_ylabel("")
    plt.xticks(rotation=45, ha="right", fontname=FONT, fontsize=10)
    plt.yticks(rotation=0, fontname=FONT, fontsize=9)

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10, pad=2)

    fig.subplots_adjust(left=0.28, right=0.90, bottom=0.14, top=0.92)

    safe_title = title_prefix.replace(" ", "_").replace(":", "").replace("/", "_")
    png_path = os.path.join(out_dir, f"{safe_title}_Heatmap{suffix}.png")
    pdf_path = os.path.join(out_dir, f"{safe_title}_Heatmap{suffix}.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {png_path}")


# ================= 6. 单模式运行 =================
def run_mode_a(ft_rel, meta, factor1, factor2):
    """
    模式A：控制 factor2，比较 factor1
    """
    mode_dir = os.path.join(OUTPUT_DIR, f"ModeA_Control_{factor2}_Compare_{factor1}")
    os.makedirs(mode_dir, exist_ok=True)

    safe_print(f"\n[START] [Mode A] Controlling {factor2}, comparing {factor1} ...")

    levels2 = natural_sort(meta[factor2].unique())

    summary_rows = []

    for lv2 in levels2:
        samples = meta.index[meta[factor2] == lv2].tolist()
        valid_samples = [s for s in samples if s in ft_rel.columns]

        if len(valid_samples) < (MIN_GROUP_SIZE * 2):
            safe_print(f"  [Skip] {factor2}={lv2}: not enough samples.")
            continue

        ft_sub = ft_rel.loc[:, valid_samples].copy()
        meta_sub = meta.loc[valid_samples].copy()

        if meta_sub[factor1].nunique() < 2:
            safe_print(f"  [Skip] {factor2}={lv2}: {factor1} has < 2 groups.")
            continue

        res_df = run_stat_analysis(ft_sub, meta_sub, compare_col=factor1)
        filt_df = filter_stat_results(res_df)

        title_prefix = f"{factor2}_{lv2}_Compare_{factor1}"
        safe_title = title_prefix.replace(" ", "_").replace(":", "").replace("/", "_")

        res_path = os.path.join(mode_dir, f"{safe_title}_Long.tsv")
        filt_path = os.path.join(mode_dir, f"{safe_title}_Filtered.tsv")
        res_df.to_csv(res_path, sep="\t", index=False)
        filt_df.to_csv(filt_path, sep="\t", index=False)

        plot_barplot(filt_df, res_df, factor1, title_prefix, mode_dir)
        plot_heatmap(filt_df, res_df, ft_sub, meta_sub, factor1, title_prefix, mode_dir)

        summary_rows.append({
            "Mode": "A",
            "Controlled_Factor": factor2,
            "Controlled_Level": lv2,
            "Compared_Factor": factor1,
            "All_Candidates": len(res_df),
            "Filtered_Candidates": len(filt_df)
        })

    if len(summary_rows) > 0:
        pd.DataFrame(summary_rows).to_csv(
            os.path.join(mode_dir, "ModeA_Summary.tsv"),
            sep="\t",
            index=False
        )


def run_mode_b(ft_rel, meta, factor1, factor2):
    """
    模式B：控制 factor1，比较 factor2
    """
    mode_dir = os.path.join(OUTPUT_DIR, f"ModeB_Control_{factor1}_Compare_{factor2}")
    os.makedirs(mode_dir, exist_ok=True)

    safe_print(f"\n[START] [Mode B] Controlling {factor1}, comparing {factor2} ...")

    levels1 = order_factor1(meta[factor1].unique())

    summary_rows = []

    for lv1 in levels1:
        samples = meta.index[meta[factor1] == lv1].tolist()
        valid_samples = [s for s in samples if s in ft_rel.columns]

        if len(valid_samples) < (MIN_GROUP_SIZE * 2):
            safe_print(f"  [Skip] {factor1}={lv1}: not enough samples.")
            continue

        ft_sub = ft_rel.loc[:, valid_samples].copy()
        meta_sub = meta.loc[valid_samples].copy()

        if meta_sub[factor2].nunique() < 2:
            safe_print(f"  [Skip] {factor1}={lv1}: {factor2} has < 2 groups.")
            continue

        res_df = run_stat_analysis(ft_sub, meta_sub, compare_col=factor2)
        filt_df = filter_stat_results(res_df)

        title_prefix = f"{factor1}_{lv1}_Compare_{factor2}"
        safe_title = title_prefix.replace(" ", "_").replace(":", "").replace("/", "_")

        res_path = os.path.join(mode_dir, f"{safe_title}_Long.tsv")
        filt_path = os.path.join(mode_dir, f"{safe_title}_Filtered.tsv")
        res_df.to_csv(res_path, sep="\t", index=False)
        filt_df.to_csv(filt_path, sep="\t", index=False)

        plot_barplot(filt_df, res_df, factor2, title_prefix, mode_dir)
        plot_heatmap(filt_df, res_df, ft_sub, meta_sub, factor2, title_prefix, mode_dir)

        summary_rows.append({
            "Mode": "B",
            "Controlled_Factor": factor1,
            "Controlled_Level": lv1,
            "Compared_Factor": factor2,
            "All_Candidates": len(res_df),
            "Filtered_Candidates": len(filt_df)
        })

    if len(summary_rows) > 0:
        pd.DataFrame(summary_rows).to_csv(
            os.path.join(mode_dir, "ModeB_Summary.tsv"),
            sep="\t",
            index=False
        )


# ================= 7. 主程序 =================
if __name__ == "__main__":
    ft_rel, meta, factor1, factor2 = load_and_clean_data()

    run_mode_a(ft_rel, meta, factor1, factor2)
    run_mode_b(ft_rel, meta, factor1, factor2)

    safe_print("\n[OK] [09_diff_viz] All bidirectional analyses finished!")