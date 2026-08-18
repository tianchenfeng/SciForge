import numpy as np
import pandas as pd

from modules.module_metabolomics.factorial import run_two_factor
from modules.module_metabolomics.multivariate import run_pca, run_plsda
from modules.module_metabolomics.preprocessing import PreprocessingConfig, preprocess_dataset
from modules.module_metabolomics.data_model import MetabolomicsDataset
from modules.module_metabolomics.statistics import bh_fdr, fold_change
from modules.module_metabolomics.univariate import run_one_factor


def test_bh_is_monotone_in_sorted_p_values():
    p = np.array([0.04, 0.001, np.nan, 0.02, 0.5])
    q = bh_fdr(p)
    ordered = np.argsort(np.nan_to_num(p, nan=np.inf))[:4]
    assert np.all(np.diff(q[ordered]) >= -1e-12)
    assert np.isnan(q[2])


def test_one_factor_detects_known_feature():
    rng = np.random.default_rng(7)
    groups = ["A"] * 6 + ["B"] * 6 + ["C"] * 6
    samples = [f"S{i}" for i in range(18)]
    metadata = pd.DataFrame({"Sample ID": samples, "Sample type": "biological", "Group": groups})
    values = rng.normal(10, .2, (5, 18)); values[0, 6:12] += 5; values[0, 12:] += 10
    abundance = pd.DataFrame(values, index=[f"F{i}" for i in range(5)], columns=samples)
    overall, contrasts = run_one_factor(abundance, metadata, "Group", "anova", "tukey", "raw")
    assert overall.set_index("Feature ID").loc["F0", "FDR"] < .05
    assert {"FC", "log2FC", "FDR", "CI low", "CI high"}.issubset(contrasts.columns)


def test_two_factor_and_simple_effects():
    rng = np.random.default_rng(11); records = []; samples = []; codes = []
    for a, a_code in [("A0", -1), ("A1", 1)]:
        for b, b_code in [("B0", -1), ("B1", 1)]:
            for replicate in range(6):
                sample = f"{a}_{b}_{replicate}"; samples.append(sample)
                records.append({"Sample ID": sample, "Sample type": "biological", "Group": a, "Diet": a, "Time": b})
                codes.append((a_code, b_code))
    metadata = pd.DataFrame(records)
    residual = rng.normal(0, .15, (6, 6)); residual = residual - residual.mean(axis=1, keepdims=True)
    data = np.tile(residual, 4)
    for index, (a, b) in enumerate(codes):
        data[:, index] += [4 * a, 4 * b, 4 * a * b, 3 * a + 3 * b, 3 * a + 3 * b + 3 * a * b, 0]
    abundance = pd.DataFrame(data, index=["A_only", "B_only", "interaction", "A_plus_B", "all", "none"], columns=samples)
    result, simple = run_two_factor(abundance, metadata, "Diet", "Time", interaction_fdr=.05)
    result = result.set_index("Feature ID")
    assert result.loc["A_only", "q_A"] < .05 and result.loc["A_only", "q_B"] > .05 and result.loc["A_only", "q_AxB"] > .05
    assert result.loc["B_only", "q_B"] < .05 and result.loc["B_only", "q_A"] > .05 and result.loc["B_only", "q_AxB"] > .05
    assert result.loc["interaction", "q_AxB"] < .05 and result.loc["interaction", "q_A"] > .05 and result.loc["interaction", "q_B"] > .05
    assert result.loc["A_plus_B", "q_A"] < .05 and result.loc["A_plus_B", "q_B"] < .05 and result.loc["A_plus_B", "q_AxB"] > .05
    assert (result.loc["all", ["q_A", "q_B", "q_AxB"]] < .05).all()
    assert (result.loc["none", ["q_A", "q_B", "q_AxB"]] > .05).all()
    assert not simple[simple["Feature ID"].eq("interaction")].empty


def test_pca_and_plsda_reproducible():
    rng = np.random.default_rng(19); samples = [f"S{i}" for i in range(16)]; groups = np.array(["A"] * 8 + ["B"] * 8)
    x = rng.normal(0, 1, (20, 16)); x[:4, groups == "B"] += 2
    abundance = pd.DataFrame(x, index=[f"F{i}" for i in range(20)], columns=samples)
    metadata = pd.DataFrame({"Sample ID": samples, "Sample type": "biological", "Group": groups})
    pca = run_pca(abundance, metadata, "pareto", 3)
    assert list(pca["scores"].columns[:4]) == ["Sample ID", "PC1", "PC2", "PC3"]
    first = run_plsda(abundance, metadata, permutations=5, random_seed=123)
    second = run_plsda(abundance, metadata, permutations=5, random_seed=123)
    assert np.allclose(first["permutation"]["Q2"], second["permutation"]["Q2"])
    assert np.isfinite(first["vip"]["VIP"]).all()


def test_random_labels_do_not_produce_excellent_plsda():
    rng = np.random.default_rng(77); samples = [f"R{i}" for i in range(24)]
    abundance = pd.DataFrame(rng.normal(size=(30, 24)), index=[f"F{i}" for i in range(30)], columns=samples)
    groups = rng.permutation(np.array(["A"] * 12 + ["B"] * 12))
    metadata = pd.DataFrame({"Sample ID": samples, "Sample type": "biological", "Group": groups})
    result = run_plsda(abundance, metadata, permutations=10, random_seed=44)
    assert result["metrics"].iloc[0]["Q2"] < .3


def test_fold_change_scales_and_explicit_preprocessing():
    raw_fc, raw_log = fold_change(np.array([2.0]), np.array([4.0]), "raw")
    log_fc, log_log = fold_change(np.array([1.0]), np.array([2.0]), "log2")
    assert raw_fc[0] == 2 and raw_log[0] == 1
    assert log_fc[0] == 2 and log_log[0] == 1
    abundance = pd.DataFrame({"S1": [0.0, 2.0], "S2": [np.nan, 4.0]}, index=["F1", "F2"])
    samples = pd.DataFrame({"Sample ID": ["S1", "S2"], "Sample type": "biological", "Group": ["A", "B"]})
    features = pd.DataFrame({"Feature ID": ["F1", "F2"]})
    dataset = MetabolomicsDataset(abundance=abundance, samples=samples, features=features, scale="raw")
    processed, _ = preprocess_dataset(dataset, PreprocessingConfig(missing_method="half_min", transformation="log2", confirmed_input_scale="raw"))
    assert processed.scale == "log2" and np.isfinite(processed.abundance.to_numpy()).all()
