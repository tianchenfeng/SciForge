import os
import re
import math
import argparse
import numpy as np
import pandas as pd
import networkx as nx
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd


# =========================================================
# 1. 路径与参数
# =========================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(description="10_export_top_genera")
parser.add_argument("--input_dir", required=True, help="01_input 文件夹路径")
parser.add_argument("--output_dir", required=True, help="02_output 文件夹路径")
parser.add_argument("--mode", required=True, choices=["single_factor", "two_factor"], help="分析模式")
args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = os.path.join(args.output_dir, "10_data_export")
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

TOP_N = 10
ALPHA = 0.05

LEVELS = {
    "Phylum": 1,
    "Genus": 5
}

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

# =========================================================
# 2. 基础工具函数
# =========================================================
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
    """
    清洗字符串，避免奇怪字符混入表头或单元格
    """
    if pd.isna(x):
        return ""
    x = str(x)
    x = x.replace("\ufeff", "")
    x = x.replace("\u200b", "")
    x = x.replace("\xa0", " ")
    x = x.strip()
    return x


def clean_header_text(x):
    """
    专门清洗表头
    """
    x = clean_text(x)
    # 去掉首尾异常中括号/空白
    x = re.sub(r"^[\[\]\)\(]+", "", x)
    x = re.sub(r"[\[\]\)\(]+$", "", x)
    x = x.strip()
    return x

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


def arcsine_sqrt_transform(percent_series):
    p = pd.to_numeric(percent_series, errors="coerce") / 100.0
    p = p.clip(lower=0, upper=1)
    return np.arcsin(np.sqrt(p))


def parse_tax_name(taxon_str, level_idx):
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


# =========================================================
# 3. 读取文件
# =========================================================
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
    factor1 = clean_header_text(meta.columns[1])   # 行因素
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
    ft = ft.apply(pd.to_numeric, errors="coerce").fillna(0)
    ft.index = ft.index.astype(str).map(clean_text)
    ft.columns = [clean_text(c) for c in ft.columns]
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


# =========================================================
# 4. 统计相关函数
# =========================================================
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
    df = df_long.dropna(subset=["Value_trans", factor1, factor2]).copy()

    if df.empty:
        return np.nan, np.nan, np.nan, False, None, "empty_data"
    if df[factor1].nunique() < 2:
        return np.nan, np.nan, np.nan, False, None, "factor1_levels_lt2"
    if df[factor2].nunique() < 2:
        return np.nan, np.nan, np.nan, False, None, "factor2_levels_lt2"
    if np.isclose(df["Value_trans"].var(ddof=1), 0, atol=1e-15):
        return np.nan, np.nan, np.nan, False, None, "zero_variance"

    try:
        formula = f"Value_trans ~ C({factor1}) + C({factor2}) + C({factor1}):C({factor2})"
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
            counts = df_long.groupby([factor1, factor2])["Value_raw"].count()
            counts = counts[counts > 0]
            if len(counts) > 0:
                avg_n = counts.mean()
                if avg_n > 0:
                    sem_model = math.sqrt(model.mse_resid / avg_n)
                    return float(sem_model), "model_based"
    except Exception:
        pass

    try:
        cell_sem = df_long.groupby([factor1, factor2])["Value_raw"].sem().dropna()
        if len(cell_sem) > 0:
            return float(cell_sem.mean()), "mean_cell_sem"
    except Exception:
        pass

    return np.nan, "NA"


def export_one_factor_table(ft_rel, top_taxa, meta, factor1, level_name):
    levels = order_row_factor(meta[factor1].tolist())
    first_col_name = clean_header_text("菌门" if level_name == "Phylum" else "菌属")
    final_rows = []
    stats_rows = []

    for taxon in top_taxa:
        taxon = clean_text(taxon)
        df_long = ft_rel.loc[[taxon]].T.join(meta[[factor1]], how="left")
        df_long = df_long.rename(columns={taxon: "Value_raw"})
        df_long["Value_raw"] = pd.to_numeric(df_long["Value_raw"], errors="coerce")
        df_long["Value_trans"] = arcsine_sqrt_transform(df_long["Value_raw"])
        df_model = df_long.dropna(subset=["Value_trans", factor1]).copy()

        model = None
        p_value = np.nan
        model_note = "ok"
        try:
            if df_model[factor1].nunique() >= 2 and not np.isclose(
                df_model["Value_trans"].var(ddof=1), 0, atol=1e-15
            ):
                model = ols(f"Value_trans ~ C({factor1})", data=df_model).fit()
                anova = sm.stats.anova_lm(model, typ=2)
                p_value = anova.loc[f"C({factor1})", "PR(>F)"]
            else:
                model_note = "insufficient_levels_or_zero_variance"
        except Exception as exc:
            model_note = f"anova_failed: {exc}"

        group_counts = df_long.groupby(factor1)["Value_raw"].count()
        avg_n = group_counts.mean() if not group_counts.empty else np.nan
        if model is not None and pd.notna(model.mse_resid) and pd.notna(avg_n) and avg_n > 0:
            sem_val = math.sqrt(model.mse_resid / avg_n)
            sem_method = "model_based_transformed_scale"
        else:
            group_sem = df_long.groupby(factor1)["Value_raw"].sem().dropna()
            sem_val = group_sem.mean() if not group_sem.empty else np.nan
            sem_method = "mean_group_sem_raw_scale" if not group_sem.empty else "NA"

        letters = {level: "" for level in levels}
        if pd.notna(p_value) and p_value < ALPHA:
            means = df_long.groupby(factor1)["Value_trans"].mean().dropna()
            groups = [level for level in levels if level in means.index]
            p_mat = run_posthoc(df_long, "Value_trans", factor1, alpha=ALPHA)
            letters.update(assign_letters_by_cliques(
                groups,
                {group: means[group] for group in groups},
                p_mat,
                alpha=ALPHA,
            ))

        for index, level in enumerate(levels):
            values = df_long.loc[df_long[factor1] == level, "Value_raw"]
            mean_value = values.mean() if not values.empty else np.nan
            final_rows.append({
                first_col_name: taxon if index == 0 else "",
                factor1: level,
                "Mean": f"{format_mean(mean_value)}{letters.get(level, '')}" if pd.notna(mean_value) else "NA",
                "SEM": format_sem(sem_val) if index == 0 else "",
                f"P_{factor1}": format_p_value(p_value) if index == 0 else "",
            })

        stats_rows.append({
            "Level": level_name,
            "Taxon": taxon,
            "Factor": factor1,
            "P_Factor": p_value,
            "SEM": sem_val,
            "SEM_Method": sem_method,
            "Model_OK": model is not None,
            "Model_Note": model_note,
            "Mean_Replicates_Per_Group": avg_n,
        })

    out_path = os.path.join(OUTPUT_DIR, f"Thesis_Table_{clean_header_text(level_name)}_WithLetters.xlsx")
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        pd.DataFrame(final_rows).to_excel(writer, sheet_name="Table", index=False)
        pd.DataFrame(stats_rows).to_excel(writer, sheet_name="Stats", index=False)
    safe_print(f"[OK] Saved: {out_path}")


# =========================================================
# 5. 主函数
# =========================================================
def export_top_taxa_tables():
    safe_print("[START] Starting export of dominant taxa tables ...")

    meta_path = find_first_existing(META_PATH_CANDIDATES)
    table_path = find_first_existing(TABLE_PATH_CANDIDATES)
    tax_path = find_first_existing(TAX_PATH_CANDIDATES)

    safe_print(f"  -> metadata: {meta_path}")
    safe_print(f"  -> feature table: {table_path}")
    safe_print(f"  -> taxonomy: {tax_path}")

    meta, factor1, factor2 = read_metadata_auto(meta_path)
    ft = read_feature_table(table_path)
    tax = read_taxonomy(tax_path)

    safe_print(f"  -> Row factor (第2列): {factor1}")
    if factor2:
        safe_print(f"  -> Column factor (第3列): {factor2}")
    else:
        safe_print("  -> Analysis mode: single factor")

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

    safe_print(f"  -> Common samples: {len(common_samples)}")
    safe_print(f"  -> Common features: {len(common_features)}")

    row_levels = order_row_factor(meta[factor1].tolist())
    col_levels = natural_sort(meta[factor2].tolist()) if factor2 else []

    safe_print(f"  -> {factor1} levels: {row_levels}")
    if factor2:
        safe_print(f"  -> {factor2} levels: {col_levels}")

    for level_name, level_idx in LEVELS.items():
        safe_print(f"\n=== Processing {level_name} ===")

        # 先从 feature 层面过滤非目标序列
        keep_features = ~tax["Taxon"].apply(is_excluded_taxon_name)

        ft_use = ft.loc[keep_features].copy()
        tax_use = tax.loc[keep_features].copy()

        if ft_use.empty:
            safe_print(f"[Warning] {level_name} 过滤后无有效 feature，跳过。")
            continue

        mapping = tax_use["Taxon"].apply(lambda x: parse_tax_name(x, level_idx))
        ft_level = ft_use.groupby(mapping).sum()

        # 再保险：聚合后的分类名继续关键词过滤
        ft_level = filter_excluded_taxa(ft_level)

        if ft_level.empty:
            safe_print(f"[Warning] {level_name} 过滤后无有效分类单元，跳过。")
            continue

        # 关键：过滤后重新计算相对丰度
        col_sums = ft_level.sum(axis=0).replace(0, np.nan)
        ft_rel = ft_level.div(col_sums, axis=1) * 100
        ft_rel = ft_rel.fillna(0)

        if ft_rel.empty:
            safe_print(f"[Warning] {level_name} 无有效数据，跳过。")
            continue

        top_taxa = (
            ft_rel.mean(axis=1)
            .sort_values(ascending=False)
            .head(TOP_N)
            .index.tolist()
        )

        if factor2 is None:
            export_one_factor_table(ft_rel, top_taxa, meta, factor1, level_name)
            continue

        final_rows = []
        stats_rows = []

        first_col_name = "菌门" if level_name == "Phylum" else "菌属"
        first_col_name = clean_header_text(first_col_name)

        for taxon in top_taxa:
            taxon = clean_text(taxon)
            safe_print(f"  -> Analyzing {taxon}")

            df_long = ft_rel.loc[[taxon]].T.join(meta[[factor1, factor2]], how="left")
            df_long = df_long.rename(columns={taxon: "Value_raw"})
            df_long["Value_raw"] = pd.to_numeric(df_long["Value_raw"], errors="coerce")
            df_long["Value_trans"] = arcsine_sqrt_transform(df_long["Value_raw"])

            p_row, p_col, p_inter, model_ok, model, model_note = fit_two_way_anova(df_long, factor1, factor2)
            sem_val, sem_method = estimate_sem(df_long, factor1, factor2, model)

            # 同一列不同处理：小写
            lower_letters = {}
            for col_lv in col_levels:
                sub = df_long[df_long[factor2] == str(col_lv)].copy()
                means = sub.groupby(factor1)["Value_trans"].mean().dropna()

                groups = [g for g in row_levels if g in means.index]
                mean_dict = {g: means[g] for g in groups}

                if len(groups) >= 2 and (
                    (pd.notna(p_row) and p_row < ALPHA) or
                    (pd.notna(p_inter) and p_inter < ALPHA)
                ):
                    p_mat = run_posthoc(sub, "Value_trans", factor1, alpha=ALPHA)
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
                means = sub.groupby(factor2)["Value_trans"].mean().dropna()

                groups = [g for g in col_levels if g in means.index]
                mean_dict = {g: means[g] for g in groups}

                if len(groups) >= 2 and (
                    (pd.notna(p_col) and p_col < ALPHA) or
                    (pd.notna(p_inter) and p_inter < ALPHA)
                ):
                    p_mat = run_posthoc(sub, "Value_trans", factor2, alpha=ALPHA)
                    res = assign_letters_by_cliques(groups, mean_dict, p_mat, alpha=ALPHA, uppercase=True)

                    if len(set(res.values())) <= 1:
                        res = {g: "" for g in groups}
                else:
                    res = {g: "" for g in groups}

                for g, l in res.items():
                    upper_letters[(clean_text(row_lv), clean_text(g))] = l

            counts_by_cell = (
                df_long.groupby([factor1, factor2])["Value_raw"]
                .count()
                .reset_index(name="N")
            )
            avg_n = counts_by_cell["N"].mean() if not counts_by_cell.empty else np.nan

            stats_rows.append({
                "Level": clean_text(level_name),
                "Taxon": taxon,
                "Row_Factor": clean_text(factor1),
                "Col_Factor": clean_text(factor2),
                "P_RowFactor": p_row,
                "P_ColFactor": p_col,
                "P_Interaction": p_inter,
                "SEM": sem_val,
                "SEM_Method": clean_text(sem_method),
                "Model_OK": model_ok,
                "Model_Note": clean_text(model_note),
                "Mean_Replicates_Per_Cell": avg_n
            })

            for i, row_lv in enumerate(row_levels):
                row_lv = clean_text(row_lv)

                row = {
                    first_col_name: taxon if i == 0 else "",
                    clean_header_text(factor1): row_lv
                }

                for col_lv in col_levels:
                    col_lv = clean_text(col_lv)

                    vals = df_long[
                        (df_long[factor1] == row_lv) &
                        (df_long[factor2] == col_lv)
                    ]["Value_raw"]

                    mean_val = vals.mean() if len(vals) > 0 else np.nan

                    low = lower_letters.get((row_lv, col_lv), "")
                    up = upper_letters.get((row_lv, col_lv), "")
                    letters = f"{low}{up}"

                    if pd.isna(mean_val):
                        cell = "NA"
                    else:
                        cell = f"{format_mean(mean_val)}{letters}"

                    row[clean_header_text(col_lv)] = cell

                if i == 0:
                    row["SEM"] = format_sem(sem_val)
                    row[f"P_{clean_header_text(factor1)}"] = format_p_value(p_row)
                    row[f"P_{clean_header_text(factor2)}"] = format_p_value(p_col)
                    row[f"P_{clean_header_text(factor1)}×{clean_header_text(factor2)}"] = format_p_value(p_inter)
                else:
                    row["SEM"] = ""
                    row[f"P_{clean_header_text(factor1)}"] = ""
                    row[f"P_{clean_header_text(factor2)}"] = ""
                    row[f"P_{clean_header_text(factor1)}×{clean_header_text(factor2)}"] = ""

                final_rows.append(row)

        df_final = pd.DataFrame(final_rows)
        df_stats = pd.DataFrame(stats_rows)

        # 再保险：强制清洗所有列名
        df_final.columns = [clean_header_text(c) for c in df_final.columns]
        df_stats.columns = [clean_header_text(c) for c in df_stats.columns]

        out_path = os.path.join(OUTPUT_DIR, f"Thesis_Table_{clean_header_text(level_name)}_WithLetters.xlsx")
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            df_final.to_excel(writer, sheet_name="Table", index=False)
            df_stats.to_excel(writer, sheet_name="Stats", index=False)

        safe_print(f"[OK] Saved: {out_path}")


# =========================================================
# 6. 运行入口
# =========================================================
if __name__ == "__main__":
    export_top_taxa_tables()
