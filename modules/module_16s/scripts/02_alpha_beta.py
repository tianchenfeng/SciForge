import os
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import argparse

from scipy.stats import entropy
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import ConvexHull
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from matplotlib.patches import Ellipse

# ================= 1. 路径与参数 =================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="02_alpha_beta_analysis")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "02_alpha_beta")
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
NPG_PALETTE = [
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
    "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"
]

ALPHA = 0.05
PCOA_POINT_SIZE = 120
PERMANOVA_N_PERM = 999


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


def format_p_value(p):
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def format_mean(v):
    if pd.isna(v):
        return "NA"
    return f"{v:.2f}"


def format_sem(v):
    if pd.isna(v):
        return "NA"
    return f"{v:.3f}"


# ================= 3. 数据读取 =================
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
    factor1 = clean_header_text(meta.columns[1])   # 行因素 / 颜色
    factor2 = clean_header_text(meta.columns[2])   # 列因素 / 形状

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

# ================= 4. 统计辅助函数 =================
def run_posthoc(sub_df, value_col, group_col, alpha=0.05):
    sub_df = sub_df.dropna(subset=[value_col, group_col]).copy()

    if sub_df.empty:
        return {}
    if sub_df[group_col].nunique() < 2:
        return {}

    counts = sub_df[group_col].value_counts()
    if (counts < 2).any():
        return {}

    if np.isclose(sub_df[value_col].var(ddof=1), 0, atol=1e-15):
        return {}

    try:
        tukey = pairwise_tukeyhsd(
            endog=sub_df[value_col].astype(float),
            groups=sub_df[group_col].astype(str),
            alpha=alpha
        )

        p_mat = {}
        for row in tukey.summary().data[1:]:
            g1, g2 = clean_text(row[0]), clean_text(row[1])
            raw_p = row[3]

            if isinstance(raw_p, str):
                raw_p = raw_p.replace("<", "").replace(">", "").replace("=", "").strip()

            p = float(raw_p)
            p_mat[(g1, g2)] = p
            p_mat[(g2, g1)] = p

        return p_mat

    except Exception as e:
        safe_print(f"[Warning] Tukey failed for {group_col}: {e}")
        return {}


def assign_letters_by_cliques(group_names, mean_dict, p_values_matrix, alpha=0.05, uppercase=False):
    group_names = [clean_text(g) for g in group_names]
    n = len(group_names)

    if n == 0:
        return {}
    if n == 1:
        return {group_names[0]: ""}

    means = [mean_dict[g] for g in group_names]

    adj = np.ones((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            g1, g2 = group_names[i], group_names[j]
            p = p_values_matrix.get((g1, g2), p_values_matrix.get((g2, g1), np.nan))

            if pd.isna(p):
                adj[i, j] = adj[j, i] = 1
            elif p < alpha:
                adj[i, j] = adj[j, i] = 0
            else:
                adj[i, j] = adj[j, i] = 1

    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j] == 1:
                G.add_edge(i, j)

    cliques = list(nx.find_cliques(G))
    if not cliques:
        return {g: "" for g in group_names}

    clique_rank = []
    for idx, clq in enumerate(cliques):
        score = np.mean([means[i] for i in clq])
        clique_rank.append((idx, score, len(clq)))
    clique_rank = sorted(clique_rank, key=lambda x: (x[1], x[2]), reverse=True)

    alphabet = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ" if uppercase else "abcdefghijklmnopqrstuvwxyz")
    if len(clique_rank) > len(alphabet):
        extra = [f"X{i}" if uppercase else f"x{i}" for i in range(1, len(clique_rank) - len(alphabet) + 1)]
        alphabet.extend(extra)

    res_letters = {i: [] for i in range(n)}
    for rank, (clq_idx, _, _) in enumerate(clique_rank):
        letter = alphabet[rank]
        for node in cliques[clq_idx]:
            res_letters[node].append(letter)

    final_res = {}
    for i, g in enumerate(group_names):
        seen = []
        for l in res_letters[i]:
            if l not in seen:
                seen.append(l)
        final_res[g] = "".join(seen)

    return final_res


def fit_two_way_anova(df_long, factor1, factor2):
    df = df_long.dropna(subset=["Value", factor1, factor2]).copy()

    if df.empty:
        return np.nan, np.nan, np.nan, False, None, "empty_data"
    if df[factor1].nunique() < 2:
        return np.nan, np.nan, np.nan, False, None, "factor1_levels_lt2"
    if df[factor2].nunique() < 2:
        return np.nan, np.nan, np.nan, False, None, "factor2_levels_lt2"
    if np.isclose(df["Value"].var(ddof=1), 0, atol=1e-15):
        return np.nan, np.nan, np.nan, False, None, "zero_variance"

    try:
        formula = f"Value ~ C({factor1}) + C({factor2}) + C({factor1}):C({factor2})"
        model = ols(formula, data=df).fit()
        anova = sm.stats.anova_lm(model, typ=2)

        idx1 = f"C({factor1})"
        idx2 = f"C({factor2})"
        idx12 = f"C({factor1}):C({factor2})"

        p1 = anova.loc[idx1, "PR(>F)"] if idx1 in anova.index else np.nan
        p2 = anova.loc[idx2, "PR(>F)"] if idx2 in anova.index else np.nan
        p12 = anova.loc[idx12, "PR(>F)"] if idx12 in anova.index else np.nan

        return p1, p2, p12, True, model, "ok"

    except Exception as e:
        safe_print(f"[Warning] ANOVA failed: {e}")
        return np.nan, np.nan, np.nan, False, None, f"anova_failed: {e}"


def estimate_sem(df_long, factor1, factor2, model):
    try:
        if model is not None and hasattr(model, "mse_resid") and pd.notna(model.mse_resid):
            counts = df_long.groupby([factor1, factor2])["Value"].count()
            counts = counts[counts > 0]
            if len(counts) > 0:
                avg_n = counts.mean()
                if avg_n > 0:
                    sem_model = math.sqrt(model.mse_resid / avg_n)
                    return float(sem_model), "model_based"
    except Exception:
        pass

    try:
        cell_sem = df_long.groupby([factor1, factor2])["Value"].sem().dropna()
        if len(cell_sem) > 0:
            return float(cell_sem.mean()), "mean_cell_sem"
    except Exception:
        pass

    return np.nan, "NA"


# ================= 5. Alpha多样性计算 =================
def calculate_alpha_metrics(counts):
    counts = np.asarray(counts, dtype=float)
    total = counts.sum()

    if total <= 0:
        return {
            "Chao1": np.nan,
            "Observed_Features": np.nan,
            "Goods_coverage": np.nan,
            "Shannon": np.nan,
            "Simpson": np.nan
        }

    p = counts / total

    shannon = entropy(p)
    simpson = 1 - np.sum(p ** 2)
    observed = np.sum(counts > 0)

    singlets = np.sum(counts == 1)
    doublets = np.sum(counts == 2)
    if doublets > 0:
        chao1 = observed + (singlets ** 2) / (2 * doublets)
    else:
        chao1 = observed + (singlets * (singlets - 1)) / 2

    goods = 1 - (singlets / total)

    return {
        "Chao1": chao1,
        "Observed_Features": observed,
        "Goods_coverage": goods,
        "Shannon": shannon,
        "Simpson": simpson
    }


def build_alpha_table(df_alpha, factor1, factor2):
    metrics = ["Chao1", "Observed_Features", "Goods_coverage", "Shannon", "Simpson"]

    row_levels = order_row_factor(df_alpha[factor1].tolist())
    col_levels = natural_sort(df_alpha[factor2].tolist())

    final_rows = []
    stats_rows = []

    for metric in metrics:
        safe_print(f"  -> Building table for {metric}")

        df_long = df_alpha[[metric, factor1, factor2]].copy()
        df_long = df_long.rename(columns={metric: "Value"})

        p_row, p_col, p_inter, model_ok, model, model_note = fit_two_way_anova(df_long, factor1, factor2)
        sem_val, sem_method = estimate_sem(df_long, factor1, factor2, model)

        # 同一列不同处理：小写
        lower_letters = {}
        for col_lv in col_levels:
            sub = df_long[df_long[factor2] == str(col_lv)].copy()
            means = sub.groupby(factor1)["Value"].mean().dropna()

            groups = [g for g in row_levels if g in means.index]
            mean_dict = {g: means[g] for g in groups}

            if len(groups) >= 2 and (
                (pd.notna(p_row) and p_row < ALPHA) or
                (pd.notna(p_inter) and p_inter < ALPHA)
            ):
                p_mat = run_posthoc(sub, "Value", factor1, alpha=ALPHA)
                res = assign_letters_by_cliques(groups, mean_dict, p_mat, alpha=ALPHA, uppercase=False)
                if len(set(res.values())) <= 1:
                    res = {g: "" for g in groups}
            else:
                res = {g: "" for g in groups}

            for g, l in res.items():
                lower_letters[(clean_text(g), clean_text(col_lv))] = l

        # 同一行不同列：大写
        upper_letters = {}
        for row_lv in row_levels:
            sub = df_long[df_long[factor1] == str(row_lv)].copy()
            means = sub.groupby(factor2)["Value"].mean().dropna()

            groups = [g for g in col_levels if g in means.index]
            mean_dict = {g: means[g] for g in groups}

            if len(groups) >= 2 and (
                (pd.notna(p_col) and p_col < ALPHA) or
                (pd.notna(p_inter) and p_inter < ALPHA)
            ):
                p_mat = run_posthoc(sub, "Value", factor2, alpha=ALPHA)
                res = assign_letters_by_cliques(groups, mean_dict, p_mat, alpha=ALPHA, uppercase=True)
                if len(set(res.values())) <= 1:
                    res = {g: "" for g in groups}
            else:
                res = {g: "" for g in groups}

            for g, l in res.items():
                upper_letters[(clean_text(row_lv), clean_text(g))] = l

        counts_by_cell = (
            df_long.groupby([factor1, factor2])["Value"]
            .count()
            .reset_index(name="N")
        )
        avg_n = counts_by_cell["N"].mean() if not counts_by_cell.empty else np.nan

        stats_rows.append({
            "Metric": metric,
            "Row_Factor": factor1,
            "Col_Factor": factor2,
            "P_RowFactor": p_row,
            "P_ColFactor": p_col,
            "P_Interaction": p_inter,
            "SEM": sem_val,
            "SEM_Method": sem_method,
            "Model_OK": model_ok,
            "Model_Note": model_note,
            "Mean_Replicates_Per_Cell": avg_n
        })

        for i, row_lv in enumerate(row_levels):
            row_lv = clean_text(row_lv)
            row = {
                "Metric": metric if i == 0 else "",
                factor1: row_lv
            }

            for col_lv in col_levels:
                col_lv = clean_text(col_lv)
                vals = df_long[
                    (df_long[factor1] == row_lv) &
                    (df_long[factor2] == col_lv)
                ]["Value"]

                mean_val = vals.mean() if len(vals) > 0 else np.nan
                low = lower_letters.get((row_lv, col_lv), "")
                up = upper_letters.get((row_lv, col_lv), "")
                letters = f"{low}{up}"

                row[col_lv] = f"{format_mean(mean_val)}{letters}" if pd.notna(mean_val) else "NA"

            if i == 0:
                row["SEM"] = format_sem(sem_val)
                row[f"P_{factor1}"] = format_p_value(p_row)
                row[f"P_{factor2}"] = format_p_value(p_col)
                row[f"P_{factor1}×{factor2}"] = format_p_value(p_inter)
            else:
                row["SEM"] = ""
                row[f"P_{factor1}"] = ""
                row[f"P_{factor2}"] = ""
                row[f"P_{factor1}×{factor2}"] = ""

            final_rows.append(row)

    df_final = pd.DataFrame(final_rows)
    df_stats = pd.DataFrame(stats_rows)

    df_final.columns = [clean_header_text(c) for c in df_final.columns]
    df_stats.columns = [clean_header_text(c) for c in df_stats.columns]

    out_path = os.path.join(OUTPUT_DIR, "Alpha_Diversity_Stats_MeanSEM.xlsx")
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df_final.to_excel(writer, sheet_name="Table", index=False)
        df_stats.to_excel(writer, sheet_name="Stats", index=False)

    safe_print(f"  [OK] Saved: {out_path}")
    return df_final, df_stats


# ================= 6. Beta多样性与PCoA =================
def compute_bray_curtis(ft):
    dist_condensed = pdist(ft.T, metric="braycurtis")
    dist_matrix = squareform(dist_condensed)
    return dist_matrix, dist_condensed


def run_pcoa_from_distance(dist_matrix):
    n = dist_matrix.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H.dot(dist_matrix ** 2).dot(H)

    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    pos_eigvals = np.maximum(eigvals, 0)
    total_pos = pos_eigvals.sum()

    if total_pos <= 0:
        raise ValueError("PCoA 正特征值和为 0，无法计算解释度。")

    explained_variance = pos_eigvals / total_pos
    coords = eigvecs[:, :2] * np.sqrt(pos_eigvals[:2])

    return coords, explained_variance


def permanova_one_factor(dist_matrix, groups, n_perm=999, seed=42):
    """
    简化版单因素 PERMANOVA
    返回 pseudo-F, p-value, R2
    """
    rng = np.random.default_rng(seed)

    groups = np.asarray(groups)
    n = len(groups)
    if n != dist_matrix.shape[0]:
        raise ValueError("groups 长度与距离矩阵维度不一致。")

    unique_groups = pd.unique(groups)
    if len(unique_groups) < 2:
        return np.nan, np.nan, np.nan

    D2 = dist_matrix ** 2
    A = -0.5 * D2
    I = np.eye(n)
    J = np.ones((n, n)) / n
    G = (I - J) @ A @ (I - J)

    def _design_matrix(g):
        levels = pd.unique(g)
        cols = []
        for lv in levels:
            cols.append((g == lv).astype(float))
        X = np.column_stack(cols)
        return X

    def _pseudo_f(g):
        X = _design_matrix(g)
        H = X @ np.linalg.pinv(X.T @ X) @ X.T
        SS_between = np.trace(H @ G)
        SS_total = np.trace(G)
        SS_within = SS_total - SS_between

        df_between = len(pd.unique(g)) - 1
        df_within = n - len(pd.unique(g))

        if df_between <= 0 or df_within <= 0 or SS_within < 0:
            return np.nan, np.nan

        ms_between = SS_between / df_between
        ms_within = SS_within / df_within
        F = ms_between / ms_within if ms_within > 0 else np.nan
        R2 = SS_between / SS_total if SS_total > 0 else np.nan
        return F, R2

    F_obs, R2_obs = _pseudo_f(groups)
    if pd.isna(F_obs):
        return np.nan, np.nan, np.nan

    perm_F = []
    for _ in range(n_perm):
        g_perm = rng.permutation(groups)
        F_perm, _ = _pseudo_f(g_perm)
        if pd.notna(F_perm):
            perm_F.append(F_perm)

    if len(perm_F) == 0:
        return F_obs, np.nan, R2_obs

    p_val = (np.sum(np.array(perm_F) >= F_obs) + 1) / (len(perm_F) + 1)
    return F_obs, p_val, R2_obs


def export_permanova_results(dist_matrix, meta, factor1, factor2):
    safe_print("  -> Running PERMANOVA ...")

    results = []

    # 因素1
    F1, p1, r21 = permanova_one_factor(dist_matrix, meta[factor1].values, n_perm=PERMANOVA_N_PERM, seed=42)
    results.append({
        "Term": factor1,
        "Pseudo_F": F1,
        "R2": r21,
        "P_value": p1
    })

    # 因素2
    F2, p2, r22 = permanova_one_factor(dist_matrix, meta[factor2].values, n_perm=PERMANOVA_N_PERM, seed=43)
    results.append({
        "Term": factor2,
        "Pseudo_F": F2,
        "R2": r22,
        "P_value": p2
    })

    # 交互组合作为附加参考
    interaction_group = meta[factor1].astype(str) + "_" + meta[factor2].astype(str)
    Fi, pi, r2i = permanova_one_factor(dist_matrix, interaction_group.values, n_perm=PERMANOVA_N_PERM, seed=44)
    results.append({
        "Term": f"{factor1}×{factor2}",
        "Pseudo_F": Fi,
        "R2": r2i,
        "P_value": pi
    })

    df_perm = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, "Beta_PERMANOVA.tsv")
    df_perm.to_csv(out_path, sep="\t", index=False)

    safe_print(f"  [OK] Saved: {out_path}")
    return df_perm


def plot_pcoa(ft, meta, factor1, factor2):
    safe_print("  -> Plotting PCoA ...")

    dist_matrix, _ = compute_bray_curtis(ft)
    coords, explained_variance = run_pcoa_from_distance(dist_matrix)

    pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
    pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

    pcoa_df = pd.DataFrame(coords, index=ft.columns, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])
    pcoa_df["Hull_Group"] = pcoa_df[factor1].astype(str) + "_" + pcoa_df[factor2].astype(str)

    unique_factor1 = order_row_factor(pcoa_df[factor1].unique())
    diet_color_map = dict(zip(unique_factor1, (NPG_PALETTE * 20)[:len(unique_factor1)]))

    unique_factor2 = natural_sort(pcoa_df[factor2].unique())
    markers = ["o", "^", "s", "D", "v", "<", ">", "P", "X", "*"]
    marker_map = dict(zip(unique_factor2, markers[:len(unique_factor2)]))

    plt.figure(figsize=(10, 8))

    plt.axhline(0, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    plt.axvline(0, color="gray", linestyle="--", linewidth=1, alpha=0.5)

    # 凸包 / 连线
    for grp in pcoa_df["Hull_Group"].unique():
        sub_pts = pcoa_df[pcoa_df["Hull_Group"] == grp]
        pts = sub_pts[["PC1", "PC2"]].values

        current_factor1 = sub_pts[factor1].iloc[0]
        color = diet_color_map[current_factor1]

        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hull_pts = pts[hull.vertices]
                plt.fill(hull_pts[:, 0], hull_pts[:, 1], color=color, alpha=0.10)
                plt.plot(hull_pts[:, 0], hull_pts[:, 1], color=color, alpha=0.4, linewidth=1, linestyle="--")
            except Exception:
                pass
        elif len(pts) == 2:
            plt.plot(pts[:, 0], pts[:, 1], color=color, alpha=0.4, linewidth=1, linestyle="--")

    sns.scatterplot(
        data=pcoa_df,
        x="PC1", y="PC2",
        hue=factor1,
        style=factor2,
        palette=diet_color_map,
        markers=marker_map,
        s=PCOA_POINT_SIZE,
        edgecolor="black",
        linewidth=0.5,
        alpha=0.9
    )

    plt.title("PCoA Analysis (Bray-Curtis)", fontname=FONT, fontsize=14)
    plt.xlabel(f"PCo1 [{pc1_var:.1f}%]", fontname=FONT, fontsize=12)
    plt.ylabel(f"PCo2 [{pc2_var:.1f}%]", fontname=FONT, fontsize=12)

    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0., frameon=False)
    plt.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "PCoA_Plot_Hull_Replicates.png")
    pdf_path = os.path.join(OUTPUT_DIR, "PCoA_Plot_Hull_Replicates.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()

    pcoa_points_path = os.path.join(OUTPUT_DIR, "PCoA_Sample_Coordinates.tsv")
    pcoa_df.to_csv(pcoa_points_path, sep="\t")

    safe_print(f"  [OK] Saved: {png_path}")
    safe_print(f"  [OK] Saved: {pcoa_points_path}")

def confidence_ellipse(x, y, ax, n_std=1.5, facecolor='none', edgecolor='black', alpha=0.18, linewidth=1.5):
    """
    给一组点画协方差椭圆
    n_std=1.5 比较适合你这种重复数不多的情况，图不会太夸张
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) < 2 or len(y) < 2:
        return None

    cov = np.cov(x, y)
    if np.linalg.det(cov) <= 0:
        return None

    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)

    ell = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=theta,
        facecolor=facecolor,
        edgecolor=edgecolor,
        alpha=alpha,
        linewidth=linewidth
    )
    ax.add_patch(ell)
    return ell


def plot_pcoa_ellipse(ft, meta, factor1, factor2):
    safe_print("  -> Plotting PCoA ellipse version ...")

    dist_matrix, _ = compute_bray_curtis(ft)
    coords, explained_variance = run_pcoa_from_distance(dist_matrix)

    pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
    pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

    pcoa_df = pd.DataFrame(coords, index=ft.columns, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])
    pcoa_df["Group"] = pcoa_df[factor1].astype(str) + "_" + pcoa_df[factor2].astype(str)

    unique_factor1 = order_row_factor(pcoa_df[factor1].unique())
    color_map = dict(zip(unique_factor1, (NPG_PALETTE * 20)[:len(unique_factor1)]))

    unique_factor2 = natural_sort(pcoa_df[factor2].unique())
    markers = ["o", "^", "s", "D", "v", "<", ">", "P", "X", "*"]
    marker_map = dict(zip(unique_factor2, markers[:len(unique_factor2)]))

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.axhline(0, color="gray", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(0, color="gray", linestyle="--", linewidth=1, alpha=0.5)

    # 先画每个组合组的椭圆
    for grp in pcoa_df["Group"].unique():
        sub = pcoa_df[pcoa_df["Group"] == grp]
        this_factor1 = sub[factor1].iloc[0]
        color = color_map[this_factor1]

        confidence_ellipse(
            sub["PC1"], sub["PC2"], ax,
            n_std=1.5,
            facecolor=color,
            edgecolor=color,
            alpha=0.15,
            linewidth=1.2
        )

    # 再画散点
    sns.scatterplot(
        data=pcoa_df,
        x="PC1", y="PC2",
        hue=factor1,
        style=factor2,
        palette=color_map,
        markers=marker_map,
        s=PCOA_POINT_SIZE,
        edgecolor="black",
        linewidth=0.5,
        alpha=0.95,
        ax=ax
    )

    ax.set_title("PCoA Analysis (Bray-Curtis)", fontname=FONT, fontsize=14)
    ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontname=FONT, fontsize=12)
    ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontname=FONT, fontsize=12)

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0., frameon=False)
    plt.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "PCoA_Plot_Ellipse.png")
    pdf_path = os.path.join(OUTPUT_DIR, "PCoA_Plot_Ellipse.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close()

    safe_print(f"  [OK] Saved: {png_path}")

def plot_pcoa_single_factor1_level_ellipse(ft, meta, factor1, factor2):
    """
    固定 factor1，每个水平单独一张图；图内展示 factor2（椭圆版）
    """
    safe_print(f"  -> Plotting single-panel ellipse PCoA by each level of {factor1} ...")

    dist_matrix, _ = compute_bray_curtis(ft)
    coords, explained_variance = run_pcoa_from_distance(dist_matrix)

    pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
    pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

    pcoa_df = pd.DataFrame(coords, index=ft.columns, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])

    levels1 = order_row_factor(pcoa_df[factor1].unique())
    levels2 = natural_sort(pcoa_df[factor2].unique())
    color_map = dict(zip(levels2, (NPG_PALETTE * 20)[:len(levels2)]))

    x_min, x_max = pcoa_df["PC1"].min(), pcoa_df["PC1"].max()
    y_min, y_max = pcoa_df["PC2"].min(), pcoa_df["PC2"].max()

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Control_{factor1}_Show_{factor2}_Ellipse")
    os.makedirs(out_dir, exist_ok=True)

    for lv1 in levels1:
        sub = pcoa_df[pcoa_df[factor1] == lv1].copy()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 6))

        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv2 in levels2:
            ss = sub[sub[factor2] == lv2]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv2],
                edgecolor="black",
                linewidth=0.5,
                label=lv2
            )

            confidence_ellipse(
                ss["PC1"], ss["PC2"], ax,
                n_std=1.5,
                facecolor=color_map[lv2],
                edgecolor=color_map[lv2],
                alpha=0.18,
                linewidth=1.2
            )

        ax.set_title(f"PCoA ({factor1} = {lv1})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_xlim(x_min * 1.1, x_max * 1.1)
        ax.set_ylim(y_min * 1.1, y_max * 1.1)

        ax.legend(frameon=False, title=factor2, bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()

        png_path = os.path.join(out_dir, f"PCoA_{factor1}_{lv1}_Show_{factor2}_Ellipse.png")
        pdf_path = os.path.join(out_dir, f"PCoA_{factor1}_{lv1}_Show_{factor2}_Ellipse.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")


def plot_pcoa_single_factor1_level_hull(ft, meta, factor1, factor2):
    """
    固定 factor1，每个水平单独一张图；图内展示 factor2（棱角 hull 版）
    """
    safe_print(f"  -> Plotting single-panel hull PCoA by each level of {factor1} ...")

    dist_matrix, _ = compute_bray_curtis(ft)
    coords, explained_variance = run_pcoa_from_distance(dist_matrix)

    pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
    pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

    pcoa_df = pd.DataFrame(coords, index=ft.columns, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])

    levels1 = order_row_factor(pcoa_df[factor1].unique())
    levels2 = natural_sort(pcoa_df[factor2].unique())
    color_map = dict(zip(levels2, (NPG_PALETTE * 20)[:len(levels2)]))

    x_min, x_max = pcoa_df["PC1"].min(), pcoa_df["PC1"].max()
    y_min, y_max = pcoa_df["PC2"].min(), pcoa_df["PC2"].max()

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Control_{factor1}_Show_{factor2}_Hull")
    os.makedirs(out_dir, exist_ok=True)

    for lv1 in levels1:
        sub = pcoa_df[pcoa_df[factor1] == lv1].copy()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 6))

        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv2 in levels2:
            ss = sub[sub[factor2] == lv2]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv2],
                edgecolor="black",
                linewidth=0.5,
                label=lv2
            )

            pts = ss[["PC1", "PC2"]].values
            if len(pts) >= 3:
                try:
                    hull = ConvexHull(pts)
                    hull_pts = pts[hull.vertices]
                    ax.fill(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv2], alpha=0.12)
                    ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv2], alpha=0.6, linewidth=1.2)
                except Exception:
                    pass
            elif len(pts) == 2:
                ax.plot(pts[:, 0], pts[:, 1], color=color_map[lv2], alpha=0.6, linewidth=1.2)

        ax.set_title(f"PCoA ({factor1} = {lv1})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_xlim(x_min * 1.1, x_max * 1.1)
        ax.set_ylim(y_min * 1.1, y_max * 1.1)

        ax.legend(frameon=False, title=factor2, bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()

        png_path = os.path.join(out_dir, f"PCoA_{factor1}_{lv1}_Show_{factor2}_Hull.png")
        pdf_path = os.path.join(out_dir, f"PCoA_{factor1}_{lv1}_Show_{factor2}_Hull.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")


def plot_pcoa_single_factor2_level_ellipse(ft, meta, factor1, factor2):
    """
    固定 factor2，每个水平单独一张图；图内展示 factor1（椭圆版）
    """
    safe_print(f"  -> Plotting single-panel ellipse PCoA by each level of {factor2} ...")

    dist_matrix, _ = compute_bray_curtis(ft)
    coords, explained_variance = run_pcoa_from_distance(dist_matrix)

    pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
    pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

    pcoa_df = pd.DataFrame(coords, index=ft.columns, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])

    levels1 = order_row_factor(pcoa_df[factor1].unique())
    levels2 = natural_sort(pcoa_df[factor2].unique())
    color_map = dict(zip(levels1, (NPG_PALETTE * 20)[:len(levels1)]))

    x_min, x_max = pcoa_df["PC1"].min(), pcoa_df["PC1"].max()
    y_min, y_max = pcoa_df["PC2"].min(), pcoa_df["PC2"].max()

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Control_{factor2}_Show_{factor1}_Ellipse")
    os.makedirs(out_dir, exist_ok=True)

    for lv2 in levels2:
        sub = pcoa_df[pcoa_df[factor2] == lv2].copy()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 6))

        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv1 in levels1:
            ss = sub[sub[factor1] == lv1]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv1],
                edgecolor="black",
                linewidth=0.5,
                label=lv1
            )

            confidence_ellipse(
                ss["PC1"], ss["PC2"], ax,
                n_std=1.5,
                facecolor=color_map[lv1],
                edgecolor=color_map[lv1],
                alpha=0.18,
                linewidth=1.2
            )

        ax.set_title(f"PCoA ({factor2} = {lv2})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_xlim(x_min * 1.1, x_max * 1.1)
        ax.set_ylim(y_min * 1.1, y_max * 1.1)

        ax.legend(frameon=False, title=factor1, bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()

        png_path = os.path.join(out_dir, f"PCoA_{factor2}_{lv2}_Show_{factor1}_Ellipse.png")
        pdf_path = os.path.join(out_dir, f"PCoA_{factor2}_{lv2}_Show_{factor1}_Ellipse.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")


def plot_pcoa_single_factor2_level_hull(ft, meta, factor1, factor2):
    """
    固定 factor2，每个水平单独一张图；图内展示 factor1（棱角 hull 版）
    """
    safe_print(f"  -> Plotting single-panel hull PCoA by each level of {factor2} ...")

    dist_matrix, _ = compute_bray_curtis(ft)
    coords, explained_variance = run_pcoa_from_distance(dist_matrix)

    pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
    pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

    pcoa_df = pd.DataFrame(coords, index=ft.columns, columns=["PC1", "PC2"]).join(meta[[factor1, factor2]])

    levels1 = order_row_factor(pcoa_df[factor1].unique())
    levels2 = natural_sort(pcoa_df[factor2].unique())
    color_map = dict(zip(levels1, (NPG_PALETTE * 20)[:len(levels1)]))

    x_min, x_max = pcoa_df["PC1"].min(), pcoa_df["PC1"].max()
    y_min, y_max = pcoa_df["PC2"].min(), pcoa_df["PC2"].max()

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Control_{factor2}_Show_{factor1}_Hull")
    os.makedirs(out_dir, exist_ok=True)

    for lv2 in levels2:
        sub = pcoa_df[pcoa_df[factor2] == lv2].copy()
        if sub.empty:
            continue

        fig, ax = plt.subplots(figsize=(7, 6))

        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv1 in levels1:
            ss = sub[sub[factor1] == lv1]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv1],
                edgecolor="black",
                linewidth=0.5,
                label=lv1
            )

            pts = ss[["PC1", "PC2"]].values
            if len(pts) >= 3:
                try:
                    hull = ConvexHull(pts)
                    hull_pts = pts[hull.vertices]
                    ax.fill(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv1], alpha=0.12)
                    ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv1], alpha=0.6, linewidth=1.2)
                except Exception:
                    pass
            elif len(pts) == 2:
                ax.plot(pts[:, 0], pts[:, 1], color=color_map[lv1], alpha=0.6, linewidth=1.2)

        ax.set_title(f"PCoA ({factor2} = {lv2})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_xlim(x_min * 1.1, x_max * 1.1)
        ax.set_ylim(y_min * 1.1, y_max * 1.1)

        ax.legend(frameon=False, title=factor1, bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()

        png_path = os.path.join(out_dir, f"PCoA_{factor2}_{lv2}_Show_{factor1}_Hull.png")
        pdf_path = os.path.join(out_dir, f"PCoA_{factor2}_{lv2}_Show_{factor1}_Hull.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")

def plot_pcoa_single_factor1_level_ellipse_recalc(ft, meta, factor1, factor2):
    """
    固定 factor1，每个水平单独重新计算 PCoA；图内展示 factor2（椭圆版）
    """
    safe_print(f"  -> Plotting recalculated ellipse PCoA by each level of {factor1} ...")

    levels1 = order_row_factor(meta[factor1].unique())
    levels2 = natural_sort(meta[factor2].unique())
    color_map = dict(zip(levels2, (NPG_PALETTE * 20)[:len(levels2)]))

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Recalc_Control_{factor1}_Show_{factor2}_Ellipse")
    os.makedirs(out_dir, exist_ok=True)

    for lv1 in levels1:
        sub_samples = meta.index[meta[factor1] == lv1].tolist()
        if len(sub_samples) < 3:
            continue

        ft_sub = ft.loc[:, sub_samples].copy()
        meta_sub = meta.loc[sub_samples].copy()

        if meta_sub[factor2].nunique() < 2:
            continue

        dist_matrix, _ = compute_bray_curtis(ft_sub)
        coords, explained_variance = run_pcoa_from_distance(dist_matrix)

        pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
        pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

        pcoa_df = pd.DataFrame(coords, index=ft_sub.columns, columns=["PC1", "PC2"]).join(meta_sub[[factor1, factor2]])

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv2 in levels2:
            ss = pcoa_df[pcoa_df[factor2] == lv2]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv2],
                edgecolor="black",
                linewidth=0.5,
                label=lv2
            )

            confidence_ellipse(
                ss["PC1"], ss["PC2"], ax,
                n_std=1.5,
                facecolor=color_map[lv2],
                edgecolor=color_map[lv2],
                alpha=0.18,
                linewidth=1.2
            )

        ax.set_title(f"Recalculated PCoA ({factor1} = {lv1})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.legend(frameon=False, title=factor2, bbox_to_anchor=(1.02, 1), loc="upper left")

        plt.tight_layout()

        png_path = os.path.join(out_dir, f"Recalc_PCoA_{factor1}_{lv1}_Show_{factor2}_Ellipse.png")
        pdf_path = os.path.join(out_dir, f"Recalc_PCoA_{factor1}_{lv1}_Show_{factor2}_Ellipse.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")


def plot_pcoa_single_factor1_level_hull_recalc(ft, meta, factor1, factor2):
    """
    固定 factor1，每个水平单独重新计算 PCoA；图内展示 factor2（hull版）
    """
    safe_print(f"  -> Plotting recalculated hull PCoA by each level of {factor1} ...")

    levels1 = order_row_factor(meta[factor1].unique())
    levels2 = natural_sort(meta[factor2].unique())
    color_map = dict(zip(levels2, (NPG_PALETTE * 20)[:len(levels2)]))

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Recalc_Control_{factor1}_Show_{factor2}_Hull")
    os.makedirs(out_dir, exist_ok=True)

    for lv1 in levels1:
        sub_samples = meta.index[meta[factor1] == lv1].tolist()
        if len(sub_samples) < 3:
            continue

        ft_sub = ft.loc[:, sub_samples].copy()
        meta_sub = meta.loc[sub_samples].copy()

        if meta_sub[factor2].nunique() < 2:
            continue

        dist_matrix, _ = compute_bray_curtis(ft_sub)
        coords, explained_variance = run_pcoa_from_distance(dist_matrix)

        pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
        pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

        pcoa_df = pd.DataFrame(coords, index=ft_sub.columns, columns=["PC1", "PC2"]).join(meta_sub[[factor1, factor2]])

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv2 in levels2:
            ss = pcoa_df[pcoa_df[factor2] == lv2]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv2],
                edgecolor="black",
                linewidth=0.5,
                label=lv2
            )

            pts = ss[["PC1", "PC2"]].values
            if len(pts) >= 3:
                try:
                    hull = ConvexHull(pts)
                    hull_pts = pts[hull.vertices]
                    ax.fill(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv2], alpha=0.12)
                    ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv2], alpha=0.6, linewidth=1.2)
                except Exception:
                    pass
            elif len(pts) == 2:
                ax.plot(pts[:, 0], pts[:, 1], color=color_map[lv2], alpha=0.6, linewidth=1.2)

        ax.set_title(f"Recalculated PCoA ({factor1} = {lv1})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.legend(frameon=False, title=factor2, bbox_to_anchor=(1.02, 1), loc="upper left")

        plt.tight_layout()

        png_path = os.path.join(out_dir, f"Recalc_PCoA_{factor1}_{lv1}_Show_{factor2}_Hull.png")
        pdf_path = os.path.join(out_dir, f"Recalc_PCoA_{factor1}_{lv1}_Show_{factor2}_Hull.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")


def plot_pcoa_single_factor2_level_ellipse_recalc(ft, meta, factor1, factor2):
    """
    固定 factor2，每个水平单独重新计算 PCoA；图内展示 factor1（椭圆版）
    """
    safe_print(f"  -> Plotting recalculated ellipse PCoA by each level of {factor2} ...")

    levels1 = order_row_factor(meta[factor1].unique())
    levels2 = natural_sort(meta[factor2].unique())
    color_map = dict(zip(levels1, (NPG_PALETTE * 20)[:len(levels1)]))

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Recalc_Control_{factor2}_Show_{factor1}_Ellipse")
    os.makedirs(out_dir, exist_ok=True)

    for lv2 in levels2:
        sub_samples = meta.index[meta[factor2] == lv2].tolist()
        if len(sub_samples) < 3:
            continue

        ft_sub = ft.loc[:, sub_samples].copy()
        meta_sub = meta.loc[sub_samples].copy()

        if meta_sub[factor1].nunique() < 2:
            continue

        dist_matrix, _ = compute_bray_curtis(ft_sub)
        coords, explained_variance = run_pcoa_from_distance(dist_matrix)

        pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
        pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

        pcoa_df = pd.DataFrame(coords, index=ft_sub.columns, columns=["PC1", "PC2"]).join(meta_sub[[factor1, factor2]])

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv1 in levels1:
            ss = pcoa_df[pcoa_df[factor1] == lv1]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv1],
                edgecolor="black",
                linewidth=0.5,
                label=lv1
            )

            confidence_ellipse(
                ss["PC1"], ss["PC2"], ax,
                n_std=1.5,
                facecolor=color_map[lv1],
                edgecolor=color_map[lv1],
                alpha=0.18,
                linewidth=1.2
            )

        ax.set_title(f"Recalculated PCoA ({factor2} = {lv2})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.legend(frameon=False, title=factor1, bbox_to_anchor=(1.02, 1), loc="upper left")

        plt.tight_layout()

        png_path = os.path.join(out_dir, f"Recalc_PCoA_{factor2}_{lv2}_Show_{factor1}_Ellipse.png")
        pdf_path = os.path.join(out_dir, f"Recalc_PCoA_{factor2}_{lv2}_Show_{factor1}_Ellipse.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")


def plot_pcoa_single_factor2_level_hull_recalc(ft, meta, factor1, factor2):
    """
    固定 factor2，每个水平单独重新计算 PCoA；图内展示 factor1（hull版）
    """
    safe_print(f"  -> Plotting recalculated hull PCoA by each level of {factor2} ...")

    levels1 = order_row_factor(meta[factor1].unique())
    levels2 = natural_sort(meta[factor2].unique())
    color_map = dict(zip(levels1, (NPG_PALETTE * 20)[:len(levels1)]))

    out_dir = os.path.join(OUTPUT_DIR, f"PCoA_Recalc_Control_{factor2}_Show_{factor1}_Hull")
    os.makedirs(out_dir, exist_ok=True)

    for lv2 in levels2:
        sub_samples = meta.index[meta[factor2] == lv2].tolist()
        if len(sub_samples) < 3:
            continue

        ft_sub = ft.loc[:, sub_samples].copy()
        meta_sub = meta.loc[sub_samples].copy()

        if meta_sub[factor1].nunique() < 2:
            continue

        dist_matrix, _ = compute_bray_curtis(ft_sub)
        coords, explained_variance = run_pcoa_from_distance(dist_matrix)

        pc1_var = explained_variance[0] * 100 if len(explained_variance) > 0 else np.nan
        pc2_var = explained_variance[1] * 100 if len(explained_variance) > 1 else np.nan

        pcoa_df = pd.DataFrame(coords, index=ft_sub.columns, columns=["PC1", "PC2"]).join(meta_sub[[factor1, factor2]])

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)

        for lv1 in levels1:
            ss = pcoa_df[pcoa_df[factor1] == lv1]
            if ss.empty:
                continue

            ax.scatter(
                ss["PC1"], ss["PC2"],
                s=90,
                color=color_map[lv1],
                edgecolor="black",
                linewidth=0.5,
                label=lv1
            )

            pts = ss[["PC1", "PC2"]].values
            if len(pts) >= 3:
                try:
                    hull = ConvexHull(pts)
                    hull_pts = pts[hull.vertices]
                    ax.fill(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv1], alpha=0.12)
                    ax.plot(hull_pts[:, 0], hull_pts[:, 1], color=color_map[lv1], alpha=0.6, linewidth=1.2)
                except Exception:
                    pass
            elif len(pts) == 2:
                ax.plot(pts[:, 0], pts[:, 1], color=color_map[lv1], alpha=0.6, linewidth=1.2)

        ax.set_title(f"Recalculated PCoA ({factor2} = {lv2})", fontsize=13, fontname=FONT)
        ax.set_xlabel(f"PCo1 [{pc1_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.set_ylabel(f"PCo2 [{pc2_var:.1f}%]", fontsize=11, fontname=FONT)
        ax.legend(frameon=False, title=factor1, bbox_to_anchor=(1.02, 1), loc="upper left")

        plt.tight_layout()

        png_path = os.path.join(out_dir, f"Recalc_PCoA_{factor2}_{lv2}_Show_{factor1}_Hull.png")
        pdf_path = os.path.join(out_dir, f"Recalc_PCoA_{factor2}_{lv2}_Show_{factor1}_Hull.pdf")
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close()

        safe_print(f"  [OK] Saved: {png_path}")

# ================= 7. 主程序 =================
def run_alpha_beta():
    safe_print("[START] Calculating Alpha/Beta Diversity ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)
    tax_path = find_first_existing(TAX_PATH_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")
    safe_print(f"  -> taxonomy: {tax_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)
    tax = read_taxonomy(tax_path)

    common_samples = [s for s in meta.index if s in ft.columns]
    if len(common_samples) == 0:
        raise ValueError("metadata 与 feature-table 没有共同样本，请检查样本名。")

    meta = meta.loc[common_samples].copy()
    ft = ft.loc[:, common_samples].copy()

    common_features = ft.index.intersection(tax.index)
    if len(common_features) == 0:
        raise ValueError("feature-table 与 taxonomy 没有共同 Feature ID。")

    ft = ft.loc[common_features].copy()
    tax = tax.loc[common_features].copy()

    # 关键：α/β 多样性也基于过滤后的有效细菌 feature 计算
    keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)
    ft = ft.loc[keep_features].copy()

    if ft.empty:
        raise ValueError("过滤线粒体/叶绿体/未分类后没有有效 feature，无法计算 α/β 多样性。")

    safe_print(f"  -> Samples: {len(common_samples)}")
    safe_print(f"  -> Row factor (第2列): {factor1}")
    safe_print(f"  -> Column factor (第3列): {factor2}")

    # 1. 计算 alpha 指标
    safe_print("  -> Calculating alpha diversity metrics ...")
    alpha_data = []
    for sample in ft.columns:
        metrics = calculate_alpha_metrics(ft[sample].values)
        row = {"SampleID": sample}
        row.update(metrics)
        alpha_data.append(row)

    df_alpha = pd.DataFrame(alpha_data).set_index("SampleID").join(meta[[factor1, factor2]])
    alpha_raw_path = os.path.join(OUTPUT_DIR, "Alpha_Diversity_Raw.tsv")
    df_alpha.to_csv(alpha_raw_path, sep="\t")
    safe_print(f"  [OK] Saved: {alpha_raw_path}")

    # 2. alpha 统计表
    build_alpha_table(df_alpha, factor1, factor2)

    # 3. beta / PCoA
    plot_pcoa(ft, meta, factor1, factor2)  # 总图：棱角 hull
    plot_pcoa_ellipse(ft, meta, factor1, factor2)  # 总图：椭圆 ellipse

    plot_pcoa_single_factor1_level_hull(ft, meta, factor1, factor2)  # 固定因素1：棱角
    plot_pcoa_single_factor1_level_ellipse(ft, meta, factor1, factor2)  # 固定因素1：椭圆

    plot_pcoa_single_factor2_level_hull(ft, meta, factor1, factor2)  # 固定因素2：棱角
    plot_pcoa_single_factor2_level_ellipse(ft, meta, factor1, factor2)  # 固定因素2：椭圆

    # 3.1 新增：子集单独重算版分图
    plot_pcoa_single_factor1_level_hull_recalc(ft, meta, factor1, factor2)
    plot_pcoa_single_factor1_level_ellipse_recalc(ft, meta, factor1, factor2)

    plot_pcoa_single_factor2_level_hull_recalc(ft, meta, factor1, factor2)
    plot_pcoa_single_factor2_level_ellipse_recalc(ft, meta, factor1, factor2)

    # 4. PERMANOVA
    dist_matrix, _ = compute_bray_curtis(ft)
    export_permanova_results(dist_matrix, meta, factor1, factor2)

    safe_print("[OK] Alpha/Beta Analysis Done.")


if __name__ == "__main__":
    run_alpha_beta()