import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="04_differential_taxa")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "04_diff")
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
TOP_N = 20
MIN_MEAN_ABUNDANCE = 0.01   # 平均相对丰度低于 0.01% 的 genus 去掉
MIN_PREVALENCE = 0.20       # 至少 20% 样本中出现
Z_CLIP = 2                  # 热图色阶截断

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

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

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
def is_excluded_taxon_name(name):
    s = clean_text(name).lower()
    return any(k in s for k in EXCLUDE_TAXA_KEYWORDS)


def filter_excluded_taxa(df):
    keep_index = [not is_excluded_taxon_name(x) for x in df.index]
    return df.loc[keep_index].copy()


def build_group_order_first_appearance(meta, factor1, factor2):
    """
    严格按 metadata 中组合首次出现顺序
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

    bad_names = {
        "", "__", "_",
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

    # 兼容 QIIME2 metadata 里的 #q2:types 行
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


def load_data_strict():
    safe_print("[Step 1] Loading data ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)
    tax_path = find_first_existing(TAX_PATH_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")
    safe_print(f"  -> taxonomy: {tax_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)
    tax = read_taxonomy(tax_path)

    # 严格按 metadata 顺序对齐样本
    valid_samples = [sid for sid in meta.index if sid in ft.columns]
    if len(valid_samples) == 0:
        raise ValueError("metadata 与 feature-table 没有共同样本，请检查样本名。")

    meta = meta.loc[valid_samples].copy()
    ft = ft.loc[:, valid_samples].copy()

    # 对齐 feature 与 taxonomy
    common_features = ft.index.intersection(tax.index)
    if len(common_features) == 0:
        raise ValueError("feature-table 与 taxonomy 没有共同 Feature ID。")

    ft = ft.loc[common_features].copy()
    tax = tax.loc[common_features].copy()

    safe_print(f"  -> Samples aligned: {len(valid_samples)}")
    safe_print(f"  -> Features aligned: {len(common_features)}")
    safe_print(f"  -> Row factor (第2列): {factor1}")
    if factor2:
        safe_print(f"  -> Col factor (第3列): {factor2}")
    else:
        safe_print("  -> Analysis mode: single factor")

    return meta, ft, tax, factor1, factor2


# ================= 4. 属水平合并 =================
def aggregate_genus_relative_abundance(ft, tax):
    safe_print("[Step 2] Aggregating to genus level ...")

    # 先从 feature 层面过滤非目标序列
    keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)

    ft_use = ft.loc[keep_features].copy()
    tax_use = tax.loc[keep_features].copy()

    if ft_use.empty:
        raise ValueError("过滤线粒体/叶绿体/未分类后没有有效 feature。")

    mapping = tax_use["Taxon"].apply(lambda x: parse_tax_name(x, 5))

    # 先聚合到 genus
    ft_genus = ft_use.groupby(mapping).sum()

    # 再保险：聚合后的 genus 名继续关键词过滤
    ft_genus = filter_excluded_taxa(ft_genus)

    if ft_genus.empty:
        raise ValueError("过滤后没有有效 genus。")

    # 关键：过滤后重新计算相对丰度，使每个样本总和为 100%
    ft_genus = ft_genus.div(ft_genus.sum(axis=0).replace(0, np.nan), axis=1) * 100
    ft_genus = ft_genus.fillna(0)

    # 平均丰度 + 出现率双重过滤
    prevalence = (ft_genus > 0).sum(axis=1) / ft_genus.shape[1]
    ft_genus = ft_genus.loc[
        (ft_genus.mean(axis=1) > MIN_MEAN_ABUNDANCE) &
        (prevalence >= MIN_PREVALENCE)
    ].copy()

    safe_print(f"  -> Retained genera after filtering: {ft_genus.shape[0]}")
    return ft_genus


# ================= 5. 组合均值表 =================
def build_group_mean_table(ft_genus, meta, factor1, factor2):
    safe_print("[Step 3] Building group mean abundance table ...")

    if factor2 is None:
        group_col = meta[factor1].astype(str)
    else:
        group_col = meta[factor1].astype(str) + "_" + meta[factor2].astype(str)

    df_trans = ft_genus.T.copy()
    df_trans["Group"] = group_col.values

    df_mean = df_trans.groupby("Group").mean().T

    ordered_groups = build_group_order_first_appearance(meta, factor1, factor2)
    df_mean = df_mean.loc[:, [g for g in ordered_groups if g in df_mean.columns]]

    return df_mean


# ================= 6. 选高变 genus 并做Z分数 =================
def prepare_variable_genera_tables(df_mean, top_n=TOP_N):
    safe_print("[Step 4] Selecting top variable genera ...")

    if df_mean.empty:
        raise ValueError("组均值表为空，无法生成热图。")

    top_n = min(top_n, df_mean.shape[0])

    # 按组间标准差排序，选前20
    variability = df_mean.std(axis=1).sort_values(ascending=False)
    top_genera = variability.head(top_n).index.tolist()

    plot_data = df_mean.loc[top_genera].copy()

    # 行 Z-score；对 std=0 做保护
    def row_zscore(x):
        sd = x.std(ddof=1)
        if pd.isna(sd) or np.isclose(sd, 0, atol=1e-15):
            return pd.Series(np.zeros(len(x)), index=x.index)
        return (x - x.mean()) / sd

    plot_z = plot_data.apply(row_zscore, axis=1)
    plot_z = plot_z.clip(lower=-Z_CLIP, upper=Z_CLIP)

    return plot_data, plot_z


# ================= 7. 绘图 =================
def plot_variable_genera_heatmap(plot_z, factor1, factor2):
    safe_print("[Step 5] Plotting variable genera heatmap ...")

    fig_height = max(7, len(plot_z) * 0.32)
    plt.figure(figsize=(10, fig_height))

    ax = sns.heatmap(
        plot_z,
        cmap="vlag",
        center=0,
        vmin=-Z_CLIP,
        vmax=Z_CLIP,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Z-score (row-scaled)", "shrink": 0.6}
    )

    ax.set_title(f"Top {plot_z.shape[0]} Variable Genera Across Groups",
                 fontsize=14)
    ax.set_xlabel("")
    ax.set_ylabel("Genus")

    plt.xticks(rotation=45, ha="right", fontname=FONT, fontsize=10)
    plt.yticks(rotation=0, fontname=FONT, fontsize=9)

    plt.tight_layout()

    out_base = os.path.join(OUTPUT_DIR, "Variable_Genera_Heatmap_Zscore")
    plt.savefig(out_base + ".pdf", dpi=300, bbox_inches="tight")
    plt.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {out_base}.png")


# ================= 8. 主程序 =================
if __name__ == "__main__":
    meta, ft, tax, factor1, factor2 = load_data_strict()

    # 1. 属水平相对丰度
    ft_genus = aggregate_genus_relative_abundance(ft, tax)
    genus_rel_path = os.path.join(OUTPUT_DIR, "Genus_RelAbundance_Sample.tsv")
    ft_genus.to_csv(genus_rel_path, sep="\t")

    # 2. 组均值
    df_mean = build_group_mean_table(ft_genus, meta, factor1, factor2)
    group_mean_path = os.path.join(OUTPUT_DIR, "Genus_RelAbundance_GroupMean.tsv")
    df_mean.to_csv(group_mean_path, sep="\t")

    # 3. 前 TOP_N 个高变 genus
    plot_data, plot_z = prepare_variable_genera_tables(df_mean, top_n=TOP_N)

    top_mean_path = os.path.join(OUTPUT_DIR, "Top_Variable_Genera_GroupMean.tsv")
    top_z_path = os.path.join(OUTPUT_DIR, "Top_Variable_Genera_Zscore.tsv")
    plot_data.to_csv(top_mean_path, sep="\t")
    plot_z.to_csv(top_z_path, sep="\t")

    # 4. 热图
    plot_variable_genera_heatmap(plot_z, factor1, factor2)

    safe_print(f"  [OK] Saved: {genus_rel_path}")
    safe_print(f"  [OK] Saved: {group_mean_path}")
    safe_print(f"  [OK] Saved: {top_mean_path}")
    safe_print(f"  [OK] Saved: {top_z_path}")

    safe_print("\n[OK] [04_diff] Done! Note: this is a variable-genera heatmap, not a formal DA test.")
