import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
except ImportError:
    sm = None
    ols = None
    pairwise_tukeyhsd = None

try:
    from modules.module_stats.letter_assign import assign_letters_by_tukey
except ImportError:
    from letter_assign import assign_letters_by_tukey


def _quote_patsy_name(column_name):
    escaped = str(column_name).replace("\\", "\\\\").replace('"', '\\"')
    return f'Q("{escaped}")'


def _rename_anova_index(anova_table, replacements):
    return anova_table.rename(index=replacements)


def _make_anova_table(rows):
    return pd.DataFrame.from_dict(
        rows,
        orient="index",
        columns=["sum_sq", "df", "F", "PR(>F)"],
    )


def _fallback_one_way_anova_table(df, factor, value_col):
    groups = [
        group[value_col].astype(float).to_numpy()
        for _, group in df.groupby(factor)
    ]
    values = df[value_col].astype(float).to_numpy()
    grand_mean = float(np.mean(values))

    ss_between = sum(len(g) * (float(np.mean(g)) - grand_mean) ** 2 for g in groups)
    ss_within = sum(float(np.sum((g - float(np.mean(g))) ** 2)) for g in groups)
    df_between = max(len(groups) - 1, 0)
    df_within = max(len(values) - len(groups), 0)
    ms_between = ss_between / df_between if df_between else np.nan
    ms_within = ss_within / df_within if df_within else np.nan
    f_value = ms_between / ms_within if ms_within and not np.isnan(ms_within) else np.nan

    return _make_anova_table({
        f"C({factor})": {
            "sum_sq": ss_between,
            "df": df_between,
            "F": f_value,
            "PR(>F)": np.nan,
        },
        "Residual": {
            "sum_sq": ss_within,
            "df": df_within,
            "F": np.nan,
            "PR(>F)": np.nan,
        },
    })


def _fallback_two_way_anova_table(df, factor1, factor2, value_col):
    y = df[value_col].astype(float)
    grand_mean = float(y.mean())
    n_total = len(df)

    factor1_means = df.groupby(factor1)[value_col].mean()
    factor2_means = df.groupby(factor2)[value_col].mean()
    cell_means = df.groupby([factor1, factor2])[value_col].mean()

    ss_factor1 = sum(
        len(df[df[factor1] == level]) * (float(mean) - grand_mean) ** 2
        for level, mean in factor1_means.items()
    )
    ss_factor2 = sum(
        len(df[df[factor2] == level]) * (float(mean) - grand_mean) ** 2
        for level, mean in factor2_means.items()
    )

    ss_cells = 0.0
    for (level1, level2), mean in cell_means.items():
        n_cell = len(df[(df[factor1] == level1) & (df[factor2] == level2)])
        ss_cells += n_cell * (float(mean) - grand_mean) ** 2

    ss_interaction = max(ss_cells - ss_factor1 - ss_factor2, 0.0)

    ss_residual = 0.0
    for (level1, level2), group in df.groupby([factor1, factor2]):
        mean = float(group[value_col].mean())
        ss_residual += float(np.sum((group[value_col].astype(float) - mean) ** 2))

    df_factor1 = max(df[factor1].nunique() - 1, 0)
    df_factor2 = max(df[factor2].nunique() - 1, 0)
    df_interaction = df_factor1 * df_factor2
    df_residual = max(n_total - df[factor1].nunique() * df[factor2].nunique(), 0)
    ms_residual = ss_residual / df_residual if df_residual else np.nan

    def f_value(ss, df_effect):
        if not df_effect or not ms_residual or np.isnan(ms_residual):
            return np.nan
        return (ss / df_effect) / ms_residual

    return _make_anova_table({
        f"C({factor1})": {
            "sum_sq": ss_factor1,
            "df": df_factor1,
            "F": f_value(ss_factor1, df_factor1),
            "PR(>F)": np.nan,
        },
        f"C({factor2})": {
            "sum_sq": ss_factor2,
            "df": df_factor2,
            "F": f_value(ss_factor2, df_factor2),
            "PR(>F)": np.nan,
        },
        f"C({factor1}):C({factor2})": {
            "sum_sq": ss_interaction,
            "df": df_interaction,
            "F": f_value(ss_interaction, df_interaction),
            "PR(>F)": np.nan,
        },
        "Residual": {
            "sum_sq": ss_residual,
            "df": df_residual,
            "F": np.nan,
            "PR(>F)": np.nan,
        },
    })


def _format_mean_sem(mean, sem, letter=""):
    if pd.isna(mean):
        return "NA"

    if pd.isna(sem):
        return f"{mean:.2f}{letter}"

    return f"{mean:.2f} +/- {sem:.2f}{letter}"


def _run_tukey_matrix(df, value_col, group_col):
    sub = df[[value_col, group_col]].dropna().copy()

    if sub[group_col].nunique() < 2:
        return {}

    counts = sub[group_col].value_counts()
    if (counts < 2).any():
        return {}

    if np.isclose(sub[value_col].var(ddof=1), 0):
        return {}

    if pairwise_tukeyhsd is None:
        return {}

    try:
        tukey = pairwise_tukeyhsd(
            endog=sub[value_col].astype(float),
            groups=sub[group_col].astype(str),
            alpha=0.05,
        )

        p_matrix = {}

        for row in tukey.summary().data[1:]:
            g1 = str(row[0])
            g2 = str(row[1])
            p = float(row[3])

            p_matrix[(g1, g2)] = p
            p_matrix[(g2, g1)] = p

        return p_matrix

    except Exception:
        return {}


def run_two_way_anova(df, factor1, factor2, value_col):
    df = df.copy()

    df[factor1] = df[factor1].astype(str)
    df[factor2] = df[factor2].astype(str)
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    df = df.dropna(subset=[factor1, factor2, value_col]).copy()

    if df.empty:
        raise ValueError("No valid data after removing missing factor/value rows.")

    q_value = _quote_patsy_name(value_col)
    q_factor1 = _quote_patsy_name(factor1)
    q_factor2 = _quote_patsy_name(factor2)

    if ols is None:
        anova_table = _fallback_two_way_anova_table(df, factor1, factor2, value_col)
    else:
        model = ols(
            f"{q_value} ~ C({q_factor1}) + C({q_factor2}) + C({q_factor1}):C({q_factor2})",
            data=df,
        ).fit()

        anova_table = sm.stats.anova_lm(model, typ=2)
        anova_table = _rename_anova_index(
            anova_table,
            {
                f"C({q_factor1})": f"C({factor1})",
                f"C({q_factor2})": f"C({factor2})",
                f"C({q_factor1}):C({q_factor2})": f"C({factor1}):C({factor2})",
            },
        )

    summary = (
        df.groupby([factor1, factor2])[value_col]
        .agg(mean="mean", sem="sem", n="count")
        .reset_index()
    )

    summary["lower_letter"] = ""

    for b_level in sorted(df[factor2].unique()):
        sub = df[df[factor2] == b_level].copy()

        group_means = sub.groupby(factor1)[value_col].mean().to_dict()
        groups = list(group_means.keys())

        p_matrix = _run_tukey_matrix(
            df=sub,
            value_col=value_col,
            group_col=factor1,
        )

        if not p_matrix:
            continue

        letters = assign_letters_by_tukey(
            groups=groups,
            mean_dict=group_means,
            p_matrix=p_matrix,
            alpha=0.05,
        )

        unique_letters = set(letters.values())

        if len(unique_letters) == 1:
            letters = {k: "" for k in letters}

        for g, letter in letters.items():
            mask = (
                (summary[factor1] == g)
                & (summary[factor2] == b_level)
            )

            summary.loc[mask, "lower_letter"] = letter

    summary["upper_letter"] = ""

    for a_level in sorted(df[factor1].unique()):
        sub = df[df[factor1] == a_level].copy()

        group_means = sub.groupby(factor2)[value_col].mean().to_dict()
        groups = list(group_means.keys())

        p_matrix = _run_tukey_matrix(
            df=sub,
            value_col=value_col,
            group_col=factor2,
        )

        if not p_matrix:
            continue

        letters = assign_letters_by_tukey(
            groups=groups,
            mean_dict=group_means,
            p_matrix=p_matrix,
            alpha=0.05,
        )

        letters = {k: v.upper() for k, v in letters.items()}
        unique_letters = set(letters.values())

        if len(unique_letters) == 1:
            letters = {k: "" for k in letters}

        for g, letter in letters.items():
            mask = (
                (summary[factor1] == a_level)
                & (summary[factor2] == g)
            )

            summary.loc[mask, "upper_letter"] = letter

    summary["final_letter"] = (
        summary["lower_letter"]
        + summary["upper_letter"]
    )

    summary["mean_sem_letter"] = summary.apply(
        lambda x: _format_mean_sem(
            x["mean"],
            x["sem"],
            x["final_letter"],
        ),
        axis=1,
    )

    return {
        "analysis_type": "two_way",
        "anova_table": anova_table,
        "summary": summary,
        "raw_df": df,
        "factor1": factor1,
        "factor2": factor2,
        "value_col": value_col,
    }


def run_one_way_anova(df, factor, value_col):
    df = df.copy()

    df[factor] = df[factor].astype(str)
    df[value_col] = pd.to_numeric(df[value_col], errors="coerce")

    df = df.dropna(subset=[factor, value_col]).copy()

    if df.empty:
        raise ValueError("No valid data after removing missing factor/value rows.")

    if df[factor].nunique() < 2:
        raise ValueError("One-way ANOVA requires at least two factor levels.")

    q_value = _quote_patsy_name(value_col)
    q_factor = _quote_patsy_name(factor)

    if ols is None:
        anova_table = _fallback_one_way_anova_table(df, factor, value_col)
    else:
        model = ols(
            f"{q_value} ~ C({q_factor})",
            data=df,
        ).fit()

        anova_table = sm.stats.anova_lm(model, typ=2)
        anova_table = _rename_anova_index(
            anova_table,
            {f"C({q_factor})": f"C({factor})"},
        )

    summary = (
        df.groupby(factor)[value_col]
        .agg(mean="mean", sem="sem", n="count")
        .reset_index()
    )

    summary["letter"] = ""

    group_means = df.groupby(factor)[value_col].mean().to_dict()
    groups = list(group_means.keys())

    p_matrix = _run_tukey_matrix(
        df=df,
        value_col=value_col,
        group_col=factor,
    )

    if p_matrix:
        letters = assign_letters_by_tukey(
            groups=groups,
            mean_dict=group_means,
            p_matrix=p_matrix,
            alpha=0.05,
        )

        if len(set(letters.values())) == 1:
            letters = {k: "" for k in letters}

        for g, letter in letters.items():
            summary.loc[summary[factor] == g, "letter"] = letter

    summary["final_letter"] = summary["letter"]

    summary["mean_sem_letter"] = summary.apply(
        lambda x: _format_mean_sem(
            x["mean"],
            x["sem"],
            x["final_letter"],
        ),
        axis=1,
    )

    return {
        "analysis_type": "one_way",
        "anova_table": anova_table,
        "summary": summary,
        "raw_df": df,
        "factor": factor,
        "value_col": value_col,
    }
