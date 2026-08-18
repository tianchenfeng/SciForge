from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from modules.module_metabolomics.statistics import bh_fdr, fold_change, holm_adjust


def _clean_groups(abundance: pd.DataFrame, metadata: pd.DataFrame, factor: str) -> tuple[list[str], list[np.ndarray]]:
    meta = metadata.set_index("Sample ID")
    bio = meta[meta["Sample type"].astype(str).str.lower().ne("qc")]
    levels = [str(v) for v in pd.unique(bio[factor]) if str(v).strip()]
    if len(levels) < 2:
        raise ValueError(f"{factor} requires at least two non-empty levels.")
    blocks = []
    for level in levels:
        ids = bio.index[bio[factor].astype(str).eq(level)].intersection(abundance.columns)
        blocks.append(abundance.loc[:, ids].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float))
    return levels, blocks


def _overall_anova(blocks: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    result = stats.f_oneway(*blocks, axis=1, nan_policy="omit")
    stacked = np.concatenate(blocks, axis=1)
    grand = np.nanmean(stacked, axis=1)
    ss_between = sum(np.sum(np.isfinite(x), axis=1) * (np.nanmean(x, axis=1) - grand) ** 2 for x in blocks)
    ss_total = np.nansum((stacked - grand[:, None]) ** 2, axis=1)
    eta2 = np.divide(ss_between, ss_total, out=np.full_like(ss_between, np.nan), where=ss_total > 0)
    df = np.full(len(stacked), len(blocks) - 1, dtype=float)
    return result.statistic, result.pvalue, df, eta2


def _overall_welch(blocks: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    means = np.column_stack([np.nanmean(x, axis=1) for x in blocks])
    variances = np.column_stack([np.nanvar(x, axis=1, ddof=1) for x in blocks])
    counts = np.column_stack([np.isfinite(x).sum(axis=1) for x in blocks]).astype(float)
    variances = np.where(variances <= 0, np.nan, variances)
    weights = counts / variances
    weight_sum = np.nansum(weights, axis=1)
    weighted_mean = np.nansum(weights * means, axis=1) / weight_sum
    k = len(blocks)
    numerator = np.nansum(weights * (means - weighted_mean[:, None]) ** 2, axis=1) / (k - 1)
    term = np.nansum(((1 - weights / weight_sum[:, None]) ** 2) / (counts - 1), axis=1)
    correction = 1 + (2 * (k - 2) / (k**2 - 1)) * term
    statistic = numerator / correction
    df2 = (k**2 - 1) / (3 * term)
    p = stats.f.sf(statistic, k - 1, df2)
    effect = np.nanvar(means, axis=1) / (np.nanvar(means, axis=1) + np.nanmean(variances, axis=1))
    return statistic, p, df2, effect


def _overall_kruskal(blocks: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    statistic = np.full(blocks[0].shape[0], np.nan)
    p = np.full_like(statistic, np.nan)
    effect = np.full_like(statistic, np.nan)
    k = len(blocks)
    for index in range(len(statistic)):
        groups = [x[index, np.isfinite(x[index])] for x in blocks]
        if any(len(x) == 0 for x in groups):
            continue
        try:
            test = stats.kruskal(*groups)
            statistic[index], p[index] = test.statistic, test.pvalue
            n = sum(map(len, groups))
            effect[index] = max((test.statistic - k + 1) / max(n - k, 1), 0)
        except ValueError:
            statistic[index], p[index], effect[index] = 0.0, 1.0, 0.0
    return statistic, p, np.full_like(statistic, k - 1), effect


def _two_group_overall(blocks: list[np.ndarray], method: str):
    a, b = blocks
    if method == "mann_whitney":
        statistic = np.full(a.shape[0], np.nan); p = statistic.copy(); effect = statistic.copy()
        for i in range(len(statistic)):
            x, y = a[i, np.isfinite(a[i])], b[i, np.isfinite(b[i])]
            if not len(x) or not len(y): continue
            test = stats.mannwhitneyu(x, y, alternative="two-sided")
            statistic[i], p[i] = test.statistic, test.pvalue
            effect[i] = 2 * test.statistic / (len(x) * len(y)) - 1
        return statistic, p, np.full_like(statistic, np.nan), effect
    test = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    n1, n2 = np.isfinite(a).sum(1), np.isfinite(b).sum(1)
    v1, v2 = np.nanvar(a, axis=1, ddof=1), np.nanvar(b, axis=1, ddof=1)
    df = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / np.maximum(n1 - 1, 1) + (v2 / n2) ** 2 / np.maximum(n2 - 1, 1))
    pooled = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / np.maximum(n1 + n2 - 2, 1))
    d = (np.nanmean(b, 1) - np.nanmean(a, 1)) / pooled
    correction = 1 - 3 / np.maximum(4 * (n1 + n2) - 9, 1)
    return test.statistic, test.pvalue, df, d * correction


def _pairwise_parametric(blocks, levels, posthoc):
    k = len(levels)
    stacked = np.concatenate(blocks, axis=1)
    means = [np.nanmean(x, 1) for x in blocks]
    counts = [np.isfinite(x).sum(1) for x in blocks]
    variances = [np.nanvar(x, 1, ddof=1) for x in blocks]
    residual_df = np.isfinite(stacked).sum(1) - k
    residual_ss = sum(np.nansum((x - np.nanmean(x, 1)[:, None]) ** 2, 1) for x in blocks)
    mse = residual_ss / np.maximum(residual_df, 1)
    pending = []
    for i, j in combinations(range(k), 2):
        difference = means[j] - means[i]
        if posthoc == "games_howell":
            se = np.sqrt(variances[i] / counts[i] + variances[j] / counts[j])
            df = se**4 / ((variances[i] / counts[i]) ** 2 / np.maximum(counts[i] - 1, 1) + (variances[j] / counts[j]) ** 2 / np.maximum(counts[j] - 1, 1))
            q = np.sqrt(2) * np.abs(difference / se)
            margin = stats.studentized_range.ppf(0.95, k, df) * se / np.sqrt(2)
            df_argument = df
        else:
            se = np.sqrt(mse / 2 * (1 / counts[i] + 1 / counts[j]))
            q = np.abs(difference) / se
            df = residual_df
            same_df = len(df) > 0 and np.all(df == df[0])
            df_argument = float(df[0]) if same_df else df
            critical = stats.studentized_range.ppf(0.95, k, df_argument)
            margin = critical * se
        pending.append((i, j, means[i], means[j], difference, q, df_argument, difference - margin, difference + margin))
    if posthoc != "games_howell" and pending and all(np.isscalar(item[6]) for item in pending):
        # SciPy's exact studentized-range integration is expensive when repeated for
        # tens of thousands of features. Evaluate the exact distribution once on a
        # dense grid and interpolate log-survival values for all balanced-design pairs.
        largest = max(float(np.nanmax(item[5])) for item in pending)
        upper = min(max(largest, 8.0), 30.0)
        grid = np.linspace(0.0, upper, int(upper * 100) + 1)
        survival = stats.studentized_range.sf(grid, k, float(pending[0][6]))
        log_survival = np.log(np.clip(survival, 1e-300, 1.0))
        probabilities = [np.exp(np.interp(np.clip(item[5], 0, upper), grid, log_survival)) for item in pending]
        for index, item in enumerate(pending):
            probabilities[index] = np.where(item[5] > upper, 0.0, probabilities[index])
    else:
        def probability(item):
            return stats.studentized_range.sf(item[5], k, item[6])
        with ThreadPoolExecutor(max_workers=min(4, len(pending))) as pool:
            probabilities = list(pool.map(probability, pending))
    return [(item[0], item[1], item[2], item[3], item[4], p, item[7], item[8]) for item, p in zip(pending, probabilities)]


def _pairwise_dunn(blocks, levels):
    outputs = []
    pairs = list(combinations(range(len(levels)), 2))
    per_pair = [[] for _ in pairs]
    for feature in range(blocks[0].shape[0]):
        groups = [x[feature, np.isfinite(x[feature])] for x in blocks]
        all_values = np.concatenate(groups)
        ranks = stats.rankdata(all_values)
        tie = stats.tiecorrect(ranks)
        n = len(all_values)
        offset = 0; rank_means = []
        for group in groups:
            rank_means.append(ranks[offset:offset + len(group)].mean()); offset += len(group)
        raw = []
        for i, j in pairs:
            variance = n * (n + 1) / 12 * tie * (1 / len(groups[i]) + 1 / len(groups[j]))
            z = (rank_means[j] - rank_means[i]) / np.sqrt(variance) if variance > 0 else 0
            raw.append(2 * stats.norm.sf(abs(z)))
        adjusted = holm_adjust(raw)
        for slot, value in zip(per_pair, adjusted): slot.append(value)
    for slot, (i, j) in zip(per_pair, pairs):
        mean_i = np.array([np.nanmean(x) for x in blocks[i]])
        mean_j = np.array([np.nanmean(x) for x in blocks[j]])
        diff = mean_j - mean_i
        outputs.append((i, j, mean_i, mean_j, diff, np.asarray(slot), np.full_like(diff, np.nan), np.full_like(diff, np.nan)))
    return outputs


def run_one_factor(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    factor: str = "Group",
    method: str = "anova",
    posthoc: str | None = None,
    scale: str = "unknown",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    levels, blocks = _clean_groups(abundance, metadata, factor)
    method = method.lower().replace("-", "_").replace(" ", "_")
    if len(levels) == 2:
        statistic, p, df, effect = _two_group_overall(blocks, "mann_whitney" if method in {"mann_whitney", "mann_whitney_u"} else "welch_t")
    elif method in {"welch", "welch_anova"}:
        statistic, p, df, effect = _overall_welch(blocks)
    elif method in {"kruskal", "kruskal_wallis"}:
        statistic, p, df, effect = _overall_kruskal(blocks)
    else:
        statistic, p, df, effect = _overall_anova(blocks)
    overall = pd.DataFrame({"Feature ID": abundance.index.astype(str), "Statistic": statistic, "df": df, "P": p, "FDR": bh_fdr(p), "effect_size": effect, "Method": method})

    if len(levels) == 2:
        posthoc = "mann_whitney" if method.startswith("mann") else "games_howell"
    posthoc = (posthoc or ("dunn" if method.startswith("kruskal") else "games_howell" if method.startswith("welch") else "tukey")).lower().replace("-", "_")
    pairwise = _pairwise_dunn(blocks, levels) if posthoc == "dunn" else _pairwise_parametric(blocks, levels, posthoc)
    frames = []
    for i, j, mean_a, mean_b, difference, pair_p, low, high in pairwise:
        fc, log2fc = fold_change(mean_a, mean_b, scale)
        frame = pd.DataFrame({
            "Contrast": f"{levels[j]} vs {levels[i]}", "Feature ID": abundance.index.astype(str),
            "Reference": levels[i], "Target": levels[j], "Mean A": mean_a, "Mean B": mean_b,
            "Difference": difference, "FC": fc, "log2FC": log2fc, "P": pair_p,
            "FDR": bh_fdr(pair_p), "CI low": low, "CI high": high, "Post-hoc": posthoc,
        })
        frame["Direction"] = np.where(frame["Difference"] > 0, "Up", np.where(frame["Difference"] < 0, "Down", "No change"))
        frames.append(frame)
    contrasts = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return overall, contrasts
