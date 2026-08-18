import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="03_taxa_composition")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "03_taxa")
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
TOP_N_PHYLA = 10
TOP_N_GENUS = 15
FIG_DPI = 300


# ================= 2. 工具函数 =================
def safe_print(msg):
    try:
        print(str(msg), flush=True)
    except UnicodeEncodeError:
        safe_msg = str(msg).encode("gbk", errors="replace").decode("gbk", errors="replace")
        print(safe_msg, flush=True)

def find_first_existing(paths):
    for p in paths:
        if os.path.exists(p):
            return p
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


def order_row_factor(values):
    vals = [clean_text(v) for v in pd.unique(pd.Series(values).astype(str))]
    priority = []
    for x in ["C", "CON", "Control", "N"]:
        if x in vals:
            priority.append(x)
    rest = [v for v in vals if v not in priority]
    return priority + natural_sort(rest)


def build_group_order(meta, factor1, factor2):
    """
    按 metadata 中组合首次出现顺序排序（严格顺序版）
    """
    if factor2 is None:
        groups = meta[factor1].astype(str)
    else:
        groups = meta[factor1].astype(str) + "_" + meta[factor2].astype(str)

    ordered_groups = []
    seen = set()

    for g in groups:
        if g not in seen:
            ordered_groups.append(g)
            seen.add(g)

    return ordered_groups


def parse_tax_name(taxon_str, level_idx):
    """
    level_idx:
      0 kingdom
      1 phylum
      2 class
      3 order
      4 family
      5 genus
      6 species
    """
    if pd.isna(taxon_str):
        return "Unclassified"

    parts = [clean_text(x) for x in str(taxon_str).split(";")]
    if level_idx >= len(parts):
        return "Unclassified"

    item = clean_text(parts[level_idx])

    if "__" in item:
        name = item.split("__", 1)[1].strip()
    else:
        name = item.strip()

    if name == "":
        return "Unclassified"

    bad_names = {
        "uncultured", "unidentified", "unknown",
        "metagenome", "ambiguous_taxa"
    }
    if name.lower() in bad_names:
        return "Unclassified"

    return clean_text(name)


# ================= 3. 数据读取 =================
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
    factor1 = clean_header_text(meta.columns[1])   # 主因素
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


# ================= 4. 丰度计算 =================
def aggregate_relative_abundance(ft, tax, level_name, level_idx):
    """
    正确流程：
    1. 先根据完整 taxonomy 字符串去掉线粒体、叶绿体、未分类等非目标序列
    2. 再聚合到目标分类层级
    3. 再按剩余有效菌群重新计算相对丰度，使每个样本总和为 100%
    """

    # 先从 feature 层面过滤，而不是聚合后再删
    keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)

    ft_use = ft.loc[keep_features].copy()
    tax_use = tax.loc[keep_features].copy()

    if ft_use.empty:
        raise ValueError(f"{level_name} 层级过滤后没有有效 feature，请检查 taxonomy。")

    mapping = tax_use["Taxon"].apply(lambda x: parse_tax_name(x, level_idx))

    ft_level = ft_use.groupby(mapping).sum()

    # 再保险：聚合后的名字里如果仍含关键词，也去掉
    keep_taxa = [not is_excluded_taxon_name(x) for x in ft_level.index]
    ft_level = ft_level.loc[keep_taxa].copy()

    if ft_level.empty:
        raise ValueError(f"{level_name} 层级过滤后没有有效分类单元。")

    # 关键：这里重新计算相对丰度，所以每个样本总和会回到 100%
    col_sums = ft_level.sum(axis=0).replace(0, np.nan)
    ft_rel = ft_level.div(col_sums, axis=1) * 100
    ft_rel = ft_rel.fillna(0)

    ft_rel.index.name = level_name
    return ft_rel


def build_group_mean_table(ft_rel, meta, factor1, factor2):
    factor_columns = [factor1] + ([factor2] if factor2 else [])
    sample_df = ft_rel.T.join(meta[factor_columns], how="left")
    if factor2 is None:
        sample_df["Group"] = sample_df[factor1].astype(str)
    else:
        sample_df["Group"] = sample_df[factor1].astype(str) + "_" + sample_df[factor2].astype(str)

    tax_cols = [c for c in sample_df.columns if c not in factor_columns + ["Group"]]
    group_mean = sample_df.groupby("Group")[tax_cols].mean().T

    ordered_groups = build_group_order(meta, factor1, factor2)
    group_mean = group_mean.loc[:, [g for g in ordered_groups if g in group_mean.columns]]

    return group_mean


def prepare_top_table_for_plot(rel_df, top_n):
    """
    输入： taxa × samples(or groups)
    输出： top taxa + Others
    """
    mean_abd = rel_df.mean(axis=1).sort_values(ascending=False)
    top_taxa = mean_abd.head(top_n).index.tolist()

    plot_df = rel_df.loc[top_taxa].copy()
    others = rel_df.drop(index=top_taxa, errors="ignore").sum(axis=0)

    plot_df.loc["Others"] = others
    plot_df = plot_df.fillna(0)

    # 再保险：Others 不允许负值
    plot_df.loc["Others"] = plot_df.loc["Others"].clip(lower=0)

    # 统一顺序：top taxa + Others
    final_order = top_taxa + ["Others"]
    plot_df = plot_df.loc[final_order]

    return plot_df, top_taxa


# ================= 5. 作图函数 =================
def get_palette(n):
    """
    给堆叠柱状图准备足够多的颜色
    """
    if n <= 20:
        return sns.color_palette("tab20", n)
    elif n <= 40:
        return sns.color_palette("tab20", 20) + sns.color_palette("tab20b", min(20, n - 20))
    else:
        return sns.color_palette("husl", n)


def plot_stacked_bar(plot_df, title, xlab, output_prefix, rotate_xticks=90):
    """
    plot_df: taxa × samples/groups
    """
    taxa = plot_df.index.tolist()
    xnames = plot_df.columns.tolist()

    colors = get_palette(len(taxa))
    color_map = dict(zip(taxa, colors))

    fig_w = max(10, len(xnames) * 0.35)
    fig_h = 6.5

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    bottom = np.zeros(len(xnames))
    x = np.arange(len(xnames))

    for taxon in taxa:
        vals = plot_df.loc[taxon].values.astype(float)
        ax.bar(
            x, vals, bottom=bottom,
            color=color_map[taxon],
            edgecolor="white",
            linewidth=0.3,
            label=taxon,
            width=0.85
        )
        bottom += vals

    ax.set_ylim(0, 100)
    ax.set_ylabel("Relative abundance (%)", fontsize=12, fontname=FONT)
    ax.set_xlabel(xlab, fontsize=12, fontname=FONT)
    ax.set_title(title, fontsize=13, fontname=FONT)

    ax.set_xticks(x)
    ax.set_xticklabels(xnames, rotation=rotate_xticks, ha="right", fontsize=8, fontname=FONT)
    yticks = ax.get_yticks()
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{int(t)}" for t in yticks], fontname=FONT)

    sns.despine(ax=ax)

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        title=""
    )

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"{output_prefix}.png"), dpi=FIG_DPI)
    plt.savefig(os.path.join(OUTPUT_DIR, f"{output_prefix}.pdf"))
    plt.close()


# ================= 6. 主流程 =================
def run_taxa_composition():
    safe_print("[START] Running taxa composition analysis ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)
    tax_path = find_first_existing(TAX_PATH_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")
    safe_print(f"  -> taxonomy: {tax_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)
    tax = read_taxonomy(tax_path)

    # 样本对齐：严格按 metadata 顺序
    common_samples = [s for s in meta.index if s in ft.columns]
    if len(common_samples) == 0:
        raise ValueError("metadata 与 feature-table 没有共同样本，请检查样本名。")

    meta = meta.loc[common_samples].copy()
    ft = ft.loc[:, common_samples].copy()

    # feature 与 taxonomy 对齐
    common_features = ft.index.intersection(tax.index)
    if len(common_features) == 0:
        raise ValueError("feature-table 与 taxonomy 没有共同 Feature ID。")

    ft = ft.loc[common_features].copy()
    tax = tax.loc[common_features].copy()

    safe_print(f"  -> Samples: {len(common_samples)}")
    safe_print(f"  -> Features: {len(common_features)}")
    safe_print(f"  -> Row factor (第2列): {factor1}")
    if factor2:
        safe_print(f"  -> Col factor (第3列): {factor2}")
    else:
        safe_print("  -> Analysis mode: single factor")

    levels = {
        "Phylum": (1, TOP_N_PHYLA),
        "Genus": (5, TOP_N_GENUS)
    }

    for level_name, (level_idx, top_n) in levels.items():
        safe_print(f"\n=== Processing {level_name} ===")

        # 1. 全样本相对丰度
        ft_rel = aggregate_relative_abundance(ft, tax, level_name, level_idx)

        sample_rel_path = os.path.join(OUTPUT_DIR, f"Sample_RelAbundance_{level_name}.tsv")
        ft_rel.to_csv(sample_rel_path, sep="\t")

        # 2. 组均值相对丰度
        group_mean = build_group_mean_table(ft_rel, meta, factor1, factor2)
        group_mean_path = os.path.join(OUTPUT_DIR, f"GroupMean_RelAbundance_{level_name}.tsv")
        group_mean.to_csv(group_mean_path, sep="\t")

        # 3. 样本图数据
        sample_plot_df, top_taxa_sample = prepare_top_table_for_plot(ft_rel, top_n)
        sample_plot_path = os.path.join(OUTPUT_DIR, f"Top_{level_name}_ForPlot_Sample.tsv")
        sample_plot_df.to_csv(sample_plot_path, sep="\t")

        # 4. 组均值图数据
        group_plot_df, top_taxa_group = prepare_top_table_for_plot(group_mean, top_n)
        group_plot_path = os.path.join(OUTPUT_DIR, f"Top_{level_name}_ForPlot_GroupMean.tsv")
        group_plot_df.to_csv(group_plot_path, sep="\t")

        # 5. 样本水平堆叠图
        plot_stacked_bar(
            sample_plot_df,
            title=f"{level_name}-level community composition (Sample level)",
            xlab="Samples",
            output_prefix=f"Composition_Sample_{level_name}",
            rotate_xticks=90
        )

        # 6. 组均值堆叠图
        plot_stacked_bar(
            group_plot_df,
            title=f"{level_name}-level community composition (Group mean)",
            xlab=f"{factor1}_{factor2}" if factor2 else factor1,
            output_prefix=f"Composition_GroupMean_{level_name}",
            rotate_xticks=45
        )

        safe_print(f"  [OK] Saved sample abundance: {sample_rel_path}")
        safe_print(f"  [OK] Saved group mean abundance: {group_mean_path}")
        safe_print(f"  [OK] Saved sample plot table: {sample_plot_path}")
        safe_print(f"  [OK] Saved group plot table: {group_plot_path}")
        safe_print(f"  [OK] Saved plots for {level_name}")

    safe_print("\n[OK] [03_taxa] All done!")


if __name__ == "__main__":
    run_taxa_composition()
