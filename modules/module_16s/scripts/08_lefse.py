import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

from scipy.stats import kruskal


# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="08_candidate_biomarkers")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "08_candidate_biomarkers")
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

ALPHA = 0.05
MIN_GROUP_SIZE = 2
MIN_MEAN_ABUNDANCE = 0.01
TOP_N_BARPLOT = 30
TOP_N_DOTPLOT = 40

LDA_LIKE_THRESHOLD = 2.0

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

LEVEL_MAP = {
    "Phylum": 1,
    "Class": 2,
    "Order": 3,
    "Family": 4,
    "Genus": 5
}

GROUP_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"
]


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


def shorten_name(name, max_len=45):
    name = clean_text(name)
    if len(name) <= max_len:
        return name
    return name[:max_len] + "..."


# ================= 3. 读取数据 =================
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
            f"metadata 至少需要 {min_columns} 列：{expected}。实际列名：{list(meta.columns)}"
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
    safe_print("[Step 1] Loading metadata / feature table / taxonomy ...")

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

    safe_print(f"  -> Samples: {len(valid_samples)}")
    safe_print(f"  -> Features: {len(common_features)}")
    safe_print(f"  -> Factor1: {factor1}")
    if factor2:
        safe_print(f"  -> Factor2: {factor2}")
    else:
        safe_print("  -> Analysis mode: single factor")

    return meta, ft, tax, factor1, factor2


# ================= 4. taxonomy 处理 =================
def parse_tax_levels(taxon_str):
    res = {k: "Unclassified" for k in LEVEL_MAP.keys()}

    if pd.isna(taxon_str):
        return res

    parts = [clean_text(x) for x in str(taxon_str).split(";")]

    for lvl_name, idx in LEVEL_MAP.items():
        if idx < len(parts):
            item = parts[idx]
            if "__" in item:
                name = item.split("__", 1)[1].strip()
            else:
                name = item.strip()

            bad = {"", "uncultured", "unidentified", "unknown", "metagenome", "ambiguous_taxa"}
            if name.lower() not in bad:
                res[lvl_name] = clean_text(name)

    return res


def build_level_abundance(ft, tax, level_name):
    # 先从 feature 层面过滤非目标序列
    keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)

    ft_use = ft.loc[keep_features].copy()
    tax_use = tax.loc[keep_features].copy()

    if ft_use.empty:
        raise ValueError(f"{level_name} 层级过滤后没有有效 feature。")

    level_series = tax_use["Taxon"].apply(lambda x: parse_tax_levels(x)[level_name])

    abd = ft_use.groupby(level_series).sum()

    # 再保险：聚合后的分类名继续关键词过滤
    abd = filter_excluded_taxa(abd)

    if abd.empty:
        return abd

    # 关键：过滤后重新计算相对丰度
    abd = abd.div(abd.sum(axis=0).replace(0, np.nan), axis=1) * 100
    abd = abd.fillna(0)

    abd = abd.loc[abd.mean(axis=1) > MIN_MEAN_ABUNDANCE].copy()
    return abd


# ================= 5. 候选标志菌筛选 =================
def compute_lda_like_score(group_means):
    diff = float(group_means.max() - group_means.min())
    score = np.log10(1 + diff) * 4
    return score


def screen_candidate_biomarkers(meta, level_abd, group_col, level_name):
    groups = natural_sort(meta[group_col].unique())
    records = []

    for taxon in level_abd.index:
        values_by_group = []
        group_means = {}
        valid_group_count = 0

        for g in groups:
            vals = level_abd.loc[taxon, meta.index[meta[group_col] == g]].values.astype(float)
            vals = vals[np.isfinite(vals)]

            if len(vals) >= MIN_GROUP_SIZE:
                values_by_group.append(vals)
                valid_group_count += 1

            group_means[g] = np.mean(vals) if len(vals) > 0 else np.nan

        if valid_group_count < 2:
            continue

        try:
            _, p = kruskal(*values_by_group)
        except Exception:
            p = np.nan

        group_means_s = pd.Series(group_means).dropna()
        if group_means_s.empty:
            continue

        dominant_group = group_means_s.idxmax()
        lda_like = compute_lda_like_score(group_means_s)

        records.append({
            "Level": level_name,
            "Taxon": taxon,
            "Kruskal_P": p,
            "Dominant_Group": dominant_group,
            "LDA_like": lda_like,
            **{f"Mean_{g}": group_means[g] for g in groups}
        })

    df = pd.DataFrame(records)
    if df.empty:
        return df, df

    df = df.sort_values(["Kruskal_P", "LDA_like"], ascending=[True, False]).reset_index(drop=True)

    df_filt = df[
        (df["Kruskal_P"] < ALPHA) &
        (df["LDA_like"] >= LDA_LIKE_THRESHOLD)
    ].copy()

    return df, df_filt


# ================= 6. 条形图 =================
def plot_candidate_biomarker_barplot(df_filt, group_col, prefix):
    if df_filt.empty:
        safe_print(f"  [Warning] {prefix}: 没有通过筛选的候选标志菌，跳过条形图。")
        return

    df_plot = df_filt.copy().sort_values("LDA_like", ascending=True).tail(TOP_N_BARPLOT)

    groups = natural_sort(df_plot["Dominant_Group"].unique())
    color_map = dict(zip(groups, (GROUP_COLORS * 20)[:len(groups)]))
    bar_colors = [color_map[g] for g in df_plot["Dominant_Group"]]

    fig_h = max(6, len(df_plot) * 0.32)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    y = np.arange(len(df_plot))
    ax.barh(y, df_plot["LDA_like"], color=bar_colors, edgecolor="black", linewidth=0.4)

    ax.set_yticks(y)
    ax.set_yticklabels([shorten_name(x, 45) for x in df_plot["Taxon"]], fontname=FONT, fontsize=9)
    ax.set_xlabel("LDA-like score", fontname=FONT, fontsize=12)
    ax.set_ylabel("")
    ax.set_title(f"Candidate biomarkers by {group_col}", fontname=FONT, fontsize=13)

    handles = []
    labels = []
    for g in groups:
        handles.append(plt.Rectangle((0, 0), 1, 1, color=color_map[g]))
        labels.append(g)
    ax.legend(handles, labels, title=group_col, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

    plt.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_LDA_Barplot.png")
    pdf_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_LDA_Barplot.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {png_path}")


# ================= 7. taxonomy path dot plot =================
def build_taxonomy_path_table(df_filt, tax):
    if df_filt.empty:
        return df_filt

    tax_path_records = []
    seen_taxa = set(df_filt["Taxon"].tolist())

    genus_tax_map = {}
    for _, row in tax.iterrows():
        levels = parse_tax_levels(row["Taxon"])
        genus = levels["Genus"]
        if genus in seen_taxa and genus not in genus_tax_map:
            genus_tax_map[genus] = levels

    for _, row in df_filt.iterrows():
        taxon = row["Taxon"]
        level = row["Level"]

        if level == "Genus" and taxon in genus_tax_map:
            levels = genus_tax_map[taxon]
        else:
            levels = {k: "Unclassified" for k in LEVEL_MAP.keys()}
            levels[level] = taxon

            matched = None
            for _, tx in tax.iterrows():
                lv = parse_tax_levels(tx["Taxon"])
                if lv.get(level, "") == taxon:
                    matched = lv
                    break
            if matched is not None:
                levels = matched
                levels[level] = taxon

        tax_path_records.append({
            "Level": level,
            "Taxon": taxon,
            "Dominant_Group": row["Dominant_Group"],
            "Kruskal_P": row["Kruskal_P"],
            "LDA_like": row["LDA_like"],
            "Phylum": levels["Phylum"],
            "Class": levels["Class"],
            "Order": levels["Order"],
            "Family": levels["Family"],
            "Genus": levels["Genus"]
        })

    out = pd.DataFrame(tax_path_records)
    out = out.sort_values(["LDA_like", "Kruskal_P"], ascending=[False, True]).reset_index(drop=True)
    return out


def plot_taxonomy_path_dotplot(df_path, group_col, prefix):
    if df_path.empty:
        safe_print(f"  [Warning] {prefix}: 没有可用于 taxonomy path dot plot 的结果，跳过。")
        return

    df_plot = df_path.head(min(TOP_N_DOTPLOT, len(df_path))).copy()

    levels = ["Phylum", "Class", "Order", "Family", "Genus"]
    x_map = {lvl: i for i, lvl in enumerate(levels)}

    y_labels = [shorten_name(x, 36) for x in df_plot["Taxon"]]
    y_pos = np.arange(len(df_plot))[::-1]

    groups = natural_sort(df_plot["Dominant_Group"].unique())
    color_map = dict(zip(groups, (GROUP_COLORS * 20)[:len(groups)]))

    lda_min = df_plot["LDA_like"].min()
    lda_max = df_plot["LDA_like"].max()

    def size_map(v):
        if lda_max == lda_min:
            return 120
        return 80 + (v - lda_min) / (lda_max - lda_min) * 220

    fig_h = max(7, len(df_plot) * 0.38)
    fig, ax = plt.subplots(figsize=(12, fig_h))

    for i, (_, row) in enumerate(df_plot.iterrows()):
        yy = y_pos[i]
        xs = [x_map[lvl] for lvl in levels]
        ax.plot(xs, [yy] * len(xs), color="#D0D0D0", linewidth=1.0, zorder=1)

    for i, (_, row) in enumerate(df_plot.iterrows()):
        yy = y_pos[i]
        color = color_map.get(row["Dominant_Group"], "#333333")
        size = size_map(row["LDA_like"])

        for lvl in levels:
            label = clean_text(row[lvl])
            if label == "" or label == "Unclassified":
                face = "#EFEFEF"
                edge = "#BEBEBE"
                alpha = 0.8
            else:
                face = color
                edge = "black"
                alpha = 0.95

            ax.scatter(
                x_map[lvl], yy,
                s=size if lvl == "Genus" else size * 0.65,
                c=face,
                edgecolor=edge,
                linewidth=0.5,
                alpha=alpha,
                zorder=3
            )

    ax.set_xticks(list(x_map.values()))
    ax.set_xticklabels(levels, fontname=FONT, fontsize=11)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontname=FONT, fontsize=9)

    ax.set_xlim(-0.5, len(levels) - 0.5)
    ax.set_ylim(-1, len(df_plot))
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"Candidate biomarkers taxonomy-path by {group_col}", fontname=FONT, fontsize=14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles_group = []
    labels_group = []
    for g in groups:
        handles_group.append(
            plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=color_map[g], markeredgecolor="black",
                       markersize=9)
        )
        labels_group.append(g)

    legend1 = ax.legend(
        handles_group, labels_group,
        title=group_col,
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )
    ax.add_artist(legend1)

    size_vals = sorted(df_plot["LDA_like"].dropna().unique())
    if len(size_vals) > 0:
        if len(size_vals) > 3:
            size_ticks = [size_vals[0], np.median(size_vals), size_vals[-1]]
        else:
            size_ticks = size_vals

        handles_size = []
        labels_size = []
        for sv in size_ticks:
            handles_size.append(
                plt.scatter([], [], s=size_map(sv), c="#999999", edgecolor="black")
            )
            labels_size.append(f"{sv:.2f}")

        ax.legend(
            handles_size, labels_size,
            title="LDA-like",
            frameon=False,
            bbox_to_anchor=(1.02, 0.55),
            loc="upper left"
        )

    plt.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_TaxonomyPath_DotPlot.png")
    pdf_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_TaxonomyPath_DotPlot.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {png_path}")


# ================= 8. 单个因素运行 =================
def run_for_group_column(meta, ft, tax, group_col, prefix):
    safe_print(f"\n========== Processing factor: {group_col} ==========")

    all_res = []
    all_filt = []

    for level_name in ["Phylum", "Class", "Order", "Family", "Genus"]:
        level_abd = build_level_abundance(ft, tax, level_name)
        df_all, df_filt = screen_candidate_biomarkers(meta, level_abd, group_col, level_name)

        if not df_all.empty:
            all_res.append(df_all)
        if not df_filt.empty:
            all_filt.append(df_filt)

    if len(all_res) == 0:
        safe_print(f"[Warning] {group_col}: 没有可分析结果。")
        return

    res_long = pd.concat(all_res, axis=0, ignore_index=True)
    filt_long = pd.concat(all_filt, axis=0, ignore_index=True) if len(all_filt) > 0 else pd.DataFrame()

    long_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_Long.tsv")
    filt_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_Filtered.tsv")
    res_long.to_csv(long_path, sep="\t", index=False)
    filt_long.to_csv(filt_path, sep="\t", index=False)

    safe_print(f"  [OK] Saved: {long_path}")
    safe_print(f"  [OK] Saved: {filt_path}")

    plot_candidate_biomarker_barplot(filt_long, group_col, prefix)

    df_path = build_taxonomy_path_table(filt_long, tax)
    path_table_path = os.path.join(OUTPUT_DIR, f"{prefix}_Candidate_Biomarkers_TaxonomyPath.tsv")
    df_path.to_csv(path_table_path, sep="\t", index=False)
    safe_print(f"  [OK] Saved: {path_table_path}")

    plot_taxonomy_path_dotplot(df_path, group_col, prefix)


# ================= 9. 主程序 =================
if __name__ == "__main__":
    safe_print("[START] Running candidate biomarker screening ...")

    meta, ft, tax, factor1, factor2 = load_data_strict()

    run_for_group_column(meta, ft, tax, factor1, prefix=factor1)
    if factor2:
        run_for_group_column(meta, ft, tax, factor2, prefix=factor2)

    safe_print("\n[OK] [08_candidate_biomarkers] Done!")
