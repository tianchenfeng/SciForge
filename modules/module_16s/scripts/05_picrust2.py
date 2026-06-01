import os
import re
import gzip
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

from sklearn.decomposition import PCA


# ================= 1. 路径与参数 =================
def get_base_dir():
    """
    源码运行：返回当前脚本目录
    exe 运行：返回 PyInstaller 解包目录
    """
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))


BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="05_picrust2_function")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
parser.add_argument("--app_base_dir", required=False, default=None, help="软件根目录/打包解包目录")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_ROOT = args.output_dir
OUTPUT_DIR = os.path.join(OUTPUT_ROOT, "05_picrust2")
ANALYSIS_MODE = args.mode
APP_BASE_DIR = args.app_base_dir if args.app_base_dir else get_base_dir()

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

NSTI_PATH_CANDIDATES = [
    os.path.join(INPUT_DIR, "picrust2", "marker_predicted_and_nsti.tsv.gz"),
    os.path.join(INPUT_DIR, "picrust2", "marker_predicted_and_nsti.tsv"),
]

FILES_CONFIG = {
    "KO": {
        "pred_candidates": [
            os.path.join(INPUT_DIR, "picrust2", "KO_predicted.tsv.gz"),
            os.path.join(INPUT_DIR, "picrust2", "KO_predicted.tsv"),
        ],
        "map_candidates": [
            os.path.join(INPUT_DIR, "picrust2", "mapping", "ko_map.tsv"),
            os.path.join(INPUT_DIR, "picrust2", "ko_map.tsv"),
            os.path.join(APP_BASE_DIR, "resources", "mapping", "ko_map.tsv"),
        ]
    },
    "EC": {
        "pred_candidates": [
            os.path.join(INPUT_DIR, "picrust2", "EC_predicted.tsv.gz"),
            os.path.join(INPUT_DIR, "picrust2", "EC_predicted.tsv"),
        ],
        "map_candidates": [
            os.path.join(INPUT_DIR, "picrust2", "mapping", "ec_map.tsv"),
            os.path.join(INPUT_DIR, "picrust2", "ec_map.tsv"),
            os.path.join(APP_BASE_DIR, "resources", "mapping", "ec_map.tsv"),
        ]
    }
}

FONT = "Arial"
NPG_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"
]

MIN_MEAN_CPM = 10
TOP_N_HEATMAP = 40
Z_CLIP = 2


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
    opener = gzip.open if file_path.endswith(".gz") else open
    with opener(file_path, "rt", encoding="utf-8-sig") as f:
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


def build_group_order_first_appearance(meta, factor1, factor2):
    groups = meta[factor1].astype(str) + "_" + meta[factor2].astype(str)
    ordered_groups = []
    seen = set()
    for g in groups:
        if g not in seen:
            ordered_groups.append(g)
            seen.add(g)
    return ordered_groups


def read_table_auto(path, index_col=0):
    sep = detect_sep(path)
    compression = "gzip" if path.endswith(".gz") else None
    df = pd.read_csv(path, sep=sep, index_col=index_col, compression=compression)
    df.index = df.index.astype(str).map(clean_text)
    df.columns = [clean_header_text(c) for c in df.columns]
    return df


# ================= 3. 读取 metadata / ASV 表 =================
def read_metadata_auto(meta_path):
    sep = detect_sep(meta_path)
    meta = pd.read_csv(meta_path, sep=sep, dtype=str, encoding="utf-8-sig")

    if meta.shape[0] > 0 and str(meta.iloc[0, 0]).startswith("#q2"):
        meta = meta.iloc[1:].copy()

    meta.columns = [clean_header_text(c) for c in meta.columns]

    if meta.shape[1] < 3:
        raise ValueError(
            f"metadata 至少需要 3 列：样本名 + 两个因素。\n当前列名：{list(meta.columns)}"
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

def load_base_data():
    safe_print("[Step 1] Loading metadata / ASV table / NSTI ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)
    tax_path = find_first_existing(TAX_PATH_CANDIDATES)
    nsti_path = find_first_existing(NSTI_PATH_CANDIDATES, allow_none=True)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")
    safe_print(f"  -> taxonomy: {tax_path}")
    safe_print(f"  -> NSTI: {nsti_path if nsti_path else 'Not found'}")
    safe_print(f"  -> APP_BASE_DIR: {APP_BASE_DIR}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)
    tax = read_taxonomy(tax_path)

    valid_samples = [sid for sid in meta.index if sid in ft.columns]
    if len(valid_samples) == 0:
        raise ValueError("metadata 与 feature-table 没有共同样本，请检查样本名。")

    meta = meta.loc[valid_samples].copy()
    ft = ft.loc[:, valid_samples].copy()

    common_features = ft.index.intersection(tax.index)
    if len(common_features) == 0:
        raise ValueError("feature-table 与 taxonomy 没有共同 Feature ID。")

    ft = ft.loc[common_features].copy()
    tax = tax.loc[common_features].copy()

    # 关键：功能预测加权前，去掉非目标 feature
    keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)
    ft = ft.loc[keep_features].copy()

    if ft.empty:
        raise ValueError("过滤线粒体/叶绿体/未分类后没有有效 feature，无法进行 PICRUSt2 功能加权。")

    nsti_df = None
    if nsti_path is not None:
        nsti_df = read_table_auto(nsti_path, index_col=0)

    safe_print(f"  -> Samples aligned: {len(valid_samples)}")
    safe_print(f"  -> Row factor (第2列): {factor1}")
    safe_print(f"  -> Col factor (第3列): {factor2}")

    return meta, ft, nsti_df, factor1, factor2


# ================= 4. 功能预测表处理 =================
def load_mapping_file(map_path):
    if map_path is None:
        return {}

    try:
        df_map = pd.read_csv(map_path, sep="\t", dtype=str)
        df_map = df_map.fillna("")
        if df_map.shape[1] < 2:
            return {}
        id_col = df_map.columns[0]
        name_col = df_map.columns[1]
        mapping = dict(zip(df_map[id_col].astype(str).map(clean_text),
                           df_map[name_col].astype(str).map(clean_text)))
        return mapping
    except Exception as e:
        safe_print(f"  [WARN] Mapping file load failed: {e}")
        return {}


def align_prediction_matrix(ft, func_pred):
    ft_asvs = set(ft.index.astype(str))
    pred_rows = set(func_pred.index.astype(str))
    pred_cols = set(func_pred.columns.astype(str))

    row_overlap = len(ft_asvs & pred_rows)
    col_overlap = len(ft_asvs & pred_cols)

    if row_overlap == 0 and col_overlap == 0:
        raise ValueError("功能预测表与 ASV 表没有交集，请检查 KO/EC predicted 文件格式。")

    if row_overlap >= col_overlap:
        pred_asv_func = func_pred.copy()
    else:
        pred_asv_func = func_pred.T.copy()

    pred_asv_func.index = pred_asv_func.index.astype(str).map(clean_text)
    pred_asv_func.columns = [clean_header_text(c) for c in pred_asv_func.columns]
    pred_asv_func = pred_asv_func.apply(pd.to_numeric, errors="coerce").fillna(0)

    common_asvs = [a for a in ft.index if a in pred_asv_func.index]
    if len(common_asvs) == 0:
        raise ValueError("ASV 对齐后仍无交集，请检查输入。")

    ft_aligned = ft.loc[common_asvs].copy()
    pred_asv_func = pred_asv_func.loc[common_asvs].copy()

    return ft_aligned, pred_asv_func


def infer_metagenome(ft_aligned, pred_asv_func):
    func_abd = ft_aligned.T.dot(pred_asv_func)
    func_abd.index.name = "SampleID"
    return func_abd


def cpm_normalize(df):
    return df.div(df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0) * 1e6


# ================= 5. NSTI =================
def plot_nsti_weighted(meta, ft, nsti_df, factor1, prefix):
    if nsti_df is None:
        safe_print("  -> NSTI file not available, skip NSTI plot.")
        return

    safe_print("  -> Calculating weighted NSTI ...")

    if "metadata_NSTI" in nsti_df.columns:
        nsti_col = "metadata_NSTI"
    else:
        nsti_col = nsti_df.columns[0]

    nsti_vals = nsti_df[nsti_col].copy()
    nsti_vals.index = nsti_vals.index.astype(str).map(clean_text)
    nsti_vals = pd.to_numeric(nsti_vals, errors="coerce")

    common_asvs = [a for a in ft.index if a in nsti_vals.index]
    if len(common_asvs) == 0:
        safe_print("  [WARN] NSTI 与 ASV 无交集，跳过 NSTI 图。")
        return

    nsti_vals = nsti_vals.loc[common_asvs]
    ft_sub = ft.loc[common_asvs].copy()

    weighted_nsti = (ft_sub.T * nsti_vals).sum(axis=1) / ft_sub.T.sum(axis=1)
    weighted_nsti = weighted_nsti.replace([np.inf, -np.inf], np.nan)

    plot_df = meta.copy()
    plot_df["Weighted_NSTI"] = weighted_nsti.reindex(plot_df.index).values

    nsti_table_path = os.path.join(OUTPUT_DIR, f"{prefix}_Weighted_NSTI.tsv")
    plot_df.to_csv(nsti_table_path, sep="\t")

    plt.figure(figsize=(6, 5))
    groups = order_row_factor(plot_df[factor1].dropna().unique())
    palette = dict(zip(groups, (NPG_COLORS * 20)[:len(groups)]))

    sns.boxplot(data=plot_df, x=factor1, y="Weighted_NSTI", hue=factor1,
                order=groups, hue_order=groups, palette=palette, width=0.55, dodge=False, legend=False)
    sns.stripplot(
        data=plot_df, x=factor1, y="Weighted_NSTI",
        order=groups, color="black", size=3, alpha=0.5
    )

    plt.axhline(0.15, color="red", linestyle="--", linewidth=1)
    plt.title(f"Weighted NSTI ({prefix})", fontname=FONT, fontsize=13)
    plt.xlabel("")
    plt.ylabel("Weighted NSTI", fontname=FONT, fontsize=12)

    sns.despine()
    plt.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, f"{prefix}_Weighted_NSTI.png")
    pdf_path = os.path.join(OUTPUT_DIR, f"{prefix}_Weighted_NSTI.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()

    safe_print(f"  [OK] Saved: {nsti_table_path}")
    safe_print(f"  [OK] Saved: {png_path}")


# ================= 6. PCA =================
def run_pca(meta, func_abd, factor1, factor2, prefix):
    safe_print(f"  -> Running PCA for {prefix} ...")

    abd_cpm = cpm_normalize(func_abd)
    abd_log = np.log1p(abd_cpm)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(abd_log.values)
    var = pca.explained_variance_ratio_

    pca_df = pd.DataFrame(coords, index=abd_log.index, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])
    pca_df["Group"] = pca_df[factor1].astype(str) + "_" + pca_df[factor2].astype(str)

    coord_path = os.path.join(OUTPUT_DIR, f"{prefix}_PCA_Coordinates.tsv")
    pca_df.to_csv(coord_path, sep="\t")

    levels1 = order_row_factor(pca_df[factor1].unique())
    levels2 = natural_sort(pca_df[factor2].unique())
    palette = dict(zip(levels1, (NPG_COLORS * 20)[:len(levels1)]))
    markers = ["o", "^", "s", "D", "v", "<", ">", "P", "X", "*"]
    marker_map = dict(zip(levels2, markers[:len(levels2)]))

    plt.figure(figsize=(8, 6))
    ax = sns.scatterplot(
        data=pca_df,
        x="PC1", y="PC2",
        hue=factor1, style=factor2,
        hue_order=levels1,
        style_order=levels2,
        palette=palette,
        markers=marker_map,
        s=90,
        edgecolor="black",
        linewidth=0.5,
        alpha=0.9
    )

    plt.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    plt.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    ax.set_title(f"PCA of inferred {prefix} profiles", fontname=FONT, fontsize=13)
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)", fontname=FONT, fontsize=12)
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)", fontname=FONT, fontsize=12)

    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    plt.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, f"{prefix}_PCA.png")
    pdf_path = os.path.join(OUTPUT_DIR, f"{prefix}_PCA.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()

    safe_print(f"  [OK] Saved: {coord_path}")
    safe_print(f"  [OK] Saved: {png_path}")


# ================= 7. Heatmap =================
def rename_function_ids(index_list, mapping):
    new_idx = []
    for fid in index_list:
        fid_clean = clean_text(fid)
        name = mapping.get(fid_clean)
        if name is None:
            clean_id = fid_clean.replace("EC:", "").replace("KO:", "")
            name = mapping.get(clean_id, fid_clean)

        if isinstance(name, str):
            if ";" in name:
                name = name.split(";")[0]
            if len(name) > 60:
                name = name[:60] + "..."
        new_idx.append(name)

    return new_idx


def plot_heatmap_all_groups(meta, func_abd, mapping, factor1, factor2, prefix):
    safe_print(f"  -> Plotting heatmap for {prefix} ...")

    abd_cpm = cpm_normalize(func_abd)
    abd_filt = abd_cpm.loc[:, abd_cpm.mean(axis=0) > MIN_MEAN_CPM].copy()

    if abd_filt.shape[1] == 0:
        safe_print(f"  [WARN] No functions passed mean CPM > {MIN_MEAN_CPM}, skip heatmap.")
        return

    meta_sub = meta.loc[abd_filt.index].copy()
    meta_sub["Group"] = meta_sub[factor1].astype(str) + "_" + meta_sub[factor2].astype(str)

    df_mean = abd_filt.join(meta_sub["Group"]).groupby("Group").mean()
    ordered_groups = build_group_order_first_appearance(meta_sub, factor1, factor2)
    df_mean = df_mean.reindex([g for g in ordered_groups if g in df_mean.index])

    n_top = min(TOP_N_HEATMAP, df_mean.shape[1])
    top_feats = df_mean.std(axis=0).sort_values(ascending=False).head(n_top).index.tolist()

    plot_data = df_mean[top_feats].T.copy()

    def row_zscore(x):
        sd = x.std(ddof=1)
        if pd.isna(sd) or np.isclose(sd, 0, atol=1e-15):
            return pd.Series(np.zeros(len(x)), index=x.index)
        return (x - x.mean()) / sd

    plot_z = plot_data.apply(row_zscore, axis=1).clip(-Z_CLIP, Z_CLIP)
    plot_z.index = rename_function_ids(plot_z.index.tolist(), mapping)

    mean_path = os.path.join(OUTPUT_DIR, f"{prefix}_Function_GroupMean.tsv")
    top_mean_path = os.path.join(OUTPUT_DIR, f"{prefix}_TopVariableFunctions_GroupMean.tsv")
    z_path = os.path.join(OUTPUT_DIR, f"{prefix}_TopVariableFunctions_Zscore.tsv")

    df_mean.to_csv(mean_path, sep="\t")
    plot_data.to_csv(top_mean_path, sep="\t")
    plot_z.to_csv(z_path, sep="\t")

    fig_h = max(6, plot_z.shape[0] * 0.35)
    fig, ax = plt.subplots(figsize=(10.8, fig_h))

    hm = sns.heatmap(
        plot_z,
        cmap="vlag",
        center=0,
        vmin=-Z_CLIP,
        vmax=Z_CLIP,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={
            "label": "Z-score",
            "shrink": 0.92,
            "pad": 0.04
        },
        ax=ax
    )

    ax.set_title(f"Top variable inferred {prefix} functions across groups",
                 fontname=FONT, fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.xticks(rotation=45, ha="right", fontname=FONT, fontsize=10)
    plt.yticks(rotation=0, fontname=FONT, fontsize=9)

    cbar = hm.collections[0].colorbar
    cbar.ax.tick_params(labelsize=10, pad=2)

    fig.subplots_adjust(left=0.33, right=0.88, bottom=0.08, top=0.95)

    png_path = os.path.join(OUTPUT_DIR, f"{prefix}_Heatmap_AllGroups.png")
    pdf_path = os.path.join(OUTPUT_DIR, f"{prefix}_Heatmap_AllGroups.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    safe_print(f"  [OK] Saved: {mean_path}")
    safe_print(f"  [OK] Saved: {top_mean_path}")
    safe_print(f"  [OK] Saved: {z_path}")
    safe_print(f"  [OK] Saved: {png_path}")


# ================= 8. 主流程 =================
def run_functional_pipeline(func_type, meta, ft, nsti_df, factor1, factor2):
    safe_print(f"\n========== Processing {func_type} ==========")

    pred_path = find_first_existing(FILES_CONFIG[func_type]["pred_candidates"], allow_none=True)
    map_path = find_first_existing(FILES_CONFIG[func_type]["map_candidates"], allow_none=True)

    if pred_path is None:
        safe_print(f"  [WARN] {func_type} prediction table not found.")
        return

    safe_print(f"  -> Prediction file: {pred_path}")
    safe_print(f"  -> Mapping file: {map_path if map_path else 'Not found'}")

    func_pred = read_table_auto(pred_path, index_col=0)
    mapping = load_mapping_file(map_path)

    ft_aligned, pred_asv_func = align_prediction_matrix(ft, func_pred)

    safe_print(f"  -> ASVs used for {func_type}: {ft_aligned.shape[0]}")
    safe_print(f"  -> Predicted functions: {pred_asv_func.shape[1]}")

    func_abd = infer_metagenome(ft_aligned, pred_asv_func)

    valid_samples = [s for s in meta.index if s in func_abd.index]
    meta_valid = meta.loc[valid_samples].copy()
    func_abd = func_abd.loc[valid_samples].copy()

    raw_path = os.path.join(OUTPUT_DIR, f"{func_type}_InferredMetagenome_Raw.tsv")
    cpm_path = os.path.join(OUTPUT_DIR, f"{func_type}_InferredMetagenome_CPM.tsv")
    func_abd.to_csv(raw_path, sep="\t")
    cpm_normalize(func_abd).to_csv(cpm_path, sep="\t")

    safe_print(f"  [OK] Saved: {raw_path}")
    safe_print(f"  [OK] Saved: {cpm_path}")

    if func_type == "KO":
        plot_nsti_weighted(meta_valid, ft_aligned, nsti_df, factor1, prefix="KO")

    run_pca(meta_valid, func_abd, factor1, factor2, prefix=func_type)
    plot_heatmap_all_groups(meta_valid, func_abd, mapping, factor1, factor2, prefix=func_type)


# ================= 9. 主程序 =================
if __name__ == "__main__":
    try:
        meta, ft, nsti_df, factor1, factor2 = load_base_data()

        run_functional_pipeline("KO", meta, ft, nsti_df, factor1, factor2)
        run_functional_pipeline("EC", meta, ft, nsti_df, factor1, factor2)

        safe_print("\n[OK] [05_picrust2] Done! Check folder '02_output/05_picrust2'")
    except Exception as e:
        safe_print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()
