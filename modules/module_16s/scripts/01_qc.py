import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="01_quality_control")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "01_qc")
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

# 稀释曲线参数
RAREFACTION_STEPS = 12
MIN_READS_FOR_RAREFACTION = 100
MIN_DEPTH_START = 100

# 绘图风格
NPG_COLORS = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"
]
FONT = "Arial"


# ================= 2. 基础工具函数 =================
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


# ================= 3. 数据读取 =================
def read_metadata_auto(meta_path):
    sep = detect_sep(meta_path)
    meta = pd.read_csv(meta_path, sep=sep, dtype=str, encoding="utf-8-sig")

    # 兼容 QIIME2 metadata 中可能存在的 #q2:types 行
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
    factor1 = clean_header_text(meta.columns[1])   # 本脚本用于分组颜色
    factor2 = None if ANALYSIS_MODE == "single_factor" else clean_header_text(meta.columns[2])

    selected_columns = [sample_col, factor1] + ([factor2] if factor2 else [])
    meta = meta[selected_columns].copy()
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

    # 去掉可能存在的注释行
    ft = ft.loc[~ft.index.astype(str).str.startswith("#"), :]
    return ft


def load_data_strict():
    safe_print("[Step 1] Loading data ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)

    # 严格按 metadata 顺序保留共有样本
    valid_samples = [sid for sid in meta.index if sid in ft.columns]
    if len(valid_samples) == 0:
        raise ValueError("metadata 与 feature-table 没有共同样本，请检查样本名。")

    meta = meta.loc[valid_samples].copy()
    ft = ft.loc[:, valid_samples].copy()

    safe_print(f"  -> Loaded {len(valid_samples)} samples.")
    safe_print(f"  -> Sample order check (top 5): {valid_samples[:5]}")
    safe_print(f"  -> Group factor (第2列): {factor1}")
    if factor2:
        safe_print(f"  -> Secondary factor (第3列): {factor2}")
    else:
        safe_print("  -> Analysis mode: single factor")

    return meta, ft, factor1, factor2


# ================= 4. QC统计表 =================
def export_qc_tables(meta, ft, group_col, second_factor):
    safe_print("\n[Step 2] Exporting QC summary tables ...")

    reads_count = ft.sum(axis=0)
    observed_features = (ft > 0).sum(axis=0)

    qc_df = meta.copy()
    qc_df["Total_Reads"] = reads_count.reindex(qc_df.index).values
    qc_df["Observed_Features"] = observed_features.reindex(qc_df.index).values

    # 每样本表
    sample_stats_path = os.path.join(OUTPUT_DIR, "QC_Sample_Stats.tsv")
    qc_df.to_csv(sample_stats_path, sep="\t")

    # 按第2列分组汇总
    group_summary = (
        qc_df.groupby(group_col)
        .agg(
            N=("Total_Reads", "size"),
            Reads_Mean=("Total_Reads", "mean"),
            Reads_SD=("Total_Reads", "std"),
            Reads_Min=("Total_Reads", "min"),
            Reads_Max=("Total_Reads", "max"),
            Features_Mean=("Observed_Features", "mean"),
            Features_SD=("Observed_Features", "std")
        )
        .reset_index()
    )

    group_summary_path = os.path.join(OUTPUT_DIR, "QC_Group_Summary.tsv")
    group_summary.to_csv(group_summary_path, sep="\t", index=False)

    # 全局统计
    overall = pd.DataFrame([{
        "Total_Samples": len(qc_df),
        "Mean_Reads": qc_df["Total_Reads"].mean(),
        "Median_Reads": qc_df["Total_Reads"].median(),
        "Min_Reads": qc_df["Total_Reads"].min(),
        "Max_Reads": qc_df["Total_Reads"].max(),
        "Mean_Observed_Features": qc_df["Observed_Features"].mean(),
        "Median_Observed_Features": qc_df["Observed_Features"].median(),
        "Group_Factor": group_col,
        "Second_Factor": second_factor
    }])

    overall_path = os.path.join(OUTPUT_DIR, "QC_Overall_Summary.tsv")
    overall.to_csv(overall_path, sep="\t", index=False)

    safe_print(f"  -> Saved: {sample_stats_path}")
    safe_print(f"  -> Saved: {group_summary_path}")
    safe_print(f"  -> Saved: {overall_path}")

    return qc_df


# ================= 5. 任务A：测序量统计图 =================
def plot_reads_summary(qc_df, group_col):
    safe_print("\n[Step 3] Plotting reads summary ...")

    plt.figure(figsize=(max(8, len(qc_df) * 0.22), 6))

    unique_groups = natural_sort(qc_df[group_col].unique())
    palette = dict(zip(unique_groups, (NPG_COLORS * 20)[:len(unique_groups)]))

    ax = sns.barplot(
        x=qc_df.index,
        y=qc_df["Total_Reads"],
        hue=qc_df[group_col],
        palette=palette,
        dodge=False,
        edgecolor="black",
        linewidth=0.5
    )

    mean_reads = qc_df["Total_Reads"].mean()
    plt.axhline(
        mean_reads,
        color="gray",
        linestyle="--",
        alpha=0.8,
        linewidth=1,
        label=f"Mean: {int(mean_reads)}"
    )

    plt.title("Sequencing Depth per Sample", fontname=FONT, fontsize=12)
    plt.ylabel("Total Reads", fontname=FONT, fontsize=12)
    plt.xlabel("")
    plt.xticks(rotation=90, fontname=FONT, fontsize=8)
    plt.yticks(fontname=FONT, fontsize=10)

    sns.despine()
    plt.legend(
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
        title=group_col
    )

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "QC_Reads_Distribution.pdf"))
    plt.savefig(os.path.join(OUTPUT_DIR, "QC_Reads_Distribution.png"), dpi=300)
    plt.close()

    safe_print("  -> Saved: QC_Reads_Distribution.png")


# ================= 6. 任务B：Observed Features统计图 =================
def plot_observed_features_summary(qc_df, group_col):
    safe_print("\n[Step 4] Plotting observed features summary ...")

    plt.figure(figsize=(max(8, len(qc_df) * 0.22), 6))

    unique_groups = natural_sort(qc_df[group_col].unique())
    palette = dict(zip(unique_groups, (NPG_COLORS * 20)[:len(unique_groups)]))

    sns.barplot(
        x=qc_df.index,
        y=qc_df["Observed_Features"],
        hue=qc_df[group_col],
        palette=palette,
        dodge=False,
        edgecolor="black",
        linewidth=0.5
    )

    mean_features = qc_df["Observed_Features"].mean()
    plt.axhline(
        mean_features,
        color="gray",
        linestyle="--",
        alpha=0.8,
        linewidth=1,
        label=f"Mean: {int(mean_features)}"
    )

    plt.title("Observed Features per Sample", fontname=FONT, fontsize=12)
    plt.ylabel("Observed Features", fontname=FONT, fontsize=12)
    plt.xlabel("")
    plt.xticks(rotation=90, fontname=FONT, fontsize=8)
    plt.yticks(fontname=FONT, fontsize=10)

    sns.despine()
    plt.legend(
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
        title=group_col
    )

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "QC_Observed_Features_Distribution.pdf"))
    plt.savefig(os.path.join(OUTPUT_DIR, "QC_Observed_Features_Distribution.png"), dpi=300)
    plt.close()

    safe_print("  -> Saved: QC_Observed_Features_Distribution.png")


# ================= 7. 任务C：稀释曲线 =================
def plot_rarefaction_curve(meta, ft, group_col, steps=RAREFACTION_STEPS):
    safe_print("\n[Step 5] Calculating rarefaction curves ...")

    counts = ft.T.copy()   # 行=样本, 列=feature
    total_reads = counts.sum(axis=1)

    valid_samples = total_reads[total_reads >= MIN_READS_FOR_RAREFACTION].index.tolist()
    if len(valid_samples) == 0:
        raise ValueError("没有样本满足稀释曲线最小 reads 阈值，请检查数据。")

    counts = counts.loc[valid_samples].copy()
    total_reads = total_reads.loc[valid_samples].copy()
    meta_sub = meta.loc[valid_samples].copy()

    # -------------------------------------------------
    # 关键修正：
    # 所有样本共用同一套 depth 网格，而不是每个样本各用各的
    # -------------------------------------------------
    global_min_reads = int(total_reads.min())
    if global_min_reads < MIN_DEPTH_START:
        start_depth = max(1, global_min_reads // 5)
    else:
        start_depth = MIN_DEPTH_START

    depth_grid = np.linspace(start_depth, global_min_reads, num=steps, dtype=int)
    depth_grid = np.unique(depth_grid)
    depth_grid = depth_grid[depth_grid > 0]

    curves_data = []

    for sample_id in counts.index:
        sample_counts = counts.loc[sample_id].astype(float)
        N = float(sample_counts.sum())
        Ni = sample_counts[sample_counts > 0].values

        if N <= 0 or len(Ni) == 0:
            continue

        for d in depth_grid:
            if d > N:
                continue

            # Hurlbert 1971 期望 richness
            # E(S_n) = S - Σ(1 - Ni/N)^n
            with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                term = np.exp(d * np.log1p(-Ni / N))
            expected_species = len(Ni) - np.sum(term)

            curves_data.append({
                "SampleID": sample_id,
                group_col: meta_sub.loc[sample_id, group_col],
                "Depth": int(d),
                "Observed_Features": float(expected_species)
            })

    if len(curves_data) == 0:
        raise ValueError("稀释曲线数据为空，请检查 feature table。")

    df_curve = pd.DataFrame(curves_data)

    # 保存原始点数据
    curve_data_path = os.path.join(OUTPUT_DIR, "QC_Rarefaction_Curve_Data.tsv")
    df_curve.to_csv(curve_data_path, sep="\t", index=False)

    # 额外保存按组汇总后的均值±SD，便于你以后直接写图或重画
    df_curve_summary = (
        df_curve.groupby([group_col, "Depth"])["Observed_Features"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "Mean", "std": "SD", "count": "N"})
    )
    summary_path = os.path.join(OUTPUT_DIR, "QC_Rarefaction_Curve_GroupSummary.tsv")
    df_curve_summary.to_csv(summary_path, sep="\t", index=False)

    # -------------------------------------------------
    # 作图：同一公共 depth 网格上按组求均值 ± SD
    # -------------------------------------------------
    plt.figure(figsize=(8, 6))

    unique_groups = natural_sort(meta_sub[group_col].unique())
    palette = dict(zip(unique_groups, (NPG_COLORS * 20)[:len(unique_groups)]))

    sns.lineplot(
        data=df_curve,
        x="Depth",
        y="Observed_Features",
        hue=group_col,
        palette=palette,
        lw=2,
        estimator="mean",
        errorbar="sd"
    )

    plt.title("Rarefaction Curve (Observed Features)", fontname=FONT, fontsize=12)
    plt.xlabel("Sequencing Depth (Reads)", fontname=FONT, fontsize=12)
    plt.ylabel("Observed Features", fontname=FONT, fontsize=12)
    plt.xticks(fontname=FONT, fontsize=10)
    plt.yticks(fontname=FONT, fontsize=10)

    sns.despine()
    plt.legend(frameon=False, title=group_col)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "QC_Rarefaction_Curve.pdf"))
    plt.savefig(os.path.join(OUTPUT_DIR, "QC_Rarefaction_Curve.png"), dpi=300)
    plt.close()

    safe_print(f"  -> Saved: {curve_data_path}")
    safe_print(f"  -> Saved: {summary_path}")
    safe_print("  -> Saved: QC_Rarefaction_Curve.png")


# ================= 8. 主程序 =================
if __name__ == "__main__":
    meta, ft, factor1, factor2 = load_data_strict()
    qc_df = export_qc_tables(meta, ft, factor1, factor2)
    plot_reads_summary(qc_df, factor1)
    plot_observed_features_summary(qc_df, factor1)
    plot_rarefaction_curve(meta, ft, factor1)
    safe_print("\n[OK] [01_qc] All done! Check folder '02_output/01_qc'")
