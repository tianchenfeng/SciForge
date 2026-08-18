import numpy as np
import pandas as pd

import modules.module_metabolomics.multivariate as multivariate
from modules.module_metabolomics.multivariate import interpret_plsda_q2, run_plsda
from modules.module_metabolomics.preprocessing import FoldPreprocessor
from modules.module_metabolomics.statistics import fold_change, fold_change_rule, summarize_significant_features


def _dataset(seed=91):
    rng = np.random.default_rng(seed)
    samples = [f"S{i}" for i in range(16)]
    groups = np.array(["A"] * 8 + ["B"] * 8)
    values = rng.normal(size=(12, 16))
    values[:3, groups == "B"] += 1.5
    abundance = pd.DataFrame(values, index=[f"F{i}" for i in range(12)], columns=samples)
    metadata = pd.DataFrame({"Sample ID": samples, "Sample type": "biological", "Group": groups})
    return abundance, metadata


def test_plsda_scaler_fit_on_outer_train_only():
    abundance, metadata = _dataset()
    result = run_plsda(abundance, metadata, permutations=0, random_seed=17)
    x = abundance[metadata["Sample ID"]].T.to_numpy()
    diagnostic = result["fold_diagnostics"][0]
    expected = x[diagnostic["train_indices"]].mean(axis=0)
    assert diagnostic["n_samples_fitted"] == len(diagnostic["train_indices"])
    assert np.allclose(diagnostic["scaler_mean"], expected)
    assert not np.allclose(diagnostic["scaler_mean"], x.mean(axis=0))


def test_plsda_inner_cv_scaling_is_fold_local():
    abundance, metadata = _dataset()
    result = run_plsda(abundance, metadata, permutations=0, random_seed=18)
    x = abundance[metadata["Sample ID"]].T.to_numpy()
    outer = result["fold_diagnostics"][0]
    outer_train = x[outer["train_indices"]]
    inner = outer["inner"][0]
    expected = outer_train[inner["train_indices"]].mean(axis=0)
    assert inner["n_samples_fitted"] == len(inner["train_indices"])
    assert np.allclose(inner["scaler_mean"], expected)
    assert not np.allclose(inner["scaler_mean"], outer_train.mean(axis=0))


def test_plsda_component_selection_does_not_use_outer_test():
    abundance, metadata = _dataset()
    first = run_plsda(abundance, metadata, permutations=0, random_seed=19)
    first_fold = first["fold_diagnostics"][0]
    changed = abundance.copy()
    test_samples = metadata.iloc[first_fold["test_indices"]]["Sample ID"]
    changed.loc[:, test_samples] += 10000
    second = run_plsda(changed, metadata, permutations=0, random_seed=19)
    assert second["fold_diagnostics"][0]["selected_components"] == first_fold["selected_components"]
    assert np.allclose(second["fold_diagnostics"][0]["scaler_mean"], first_fold["scaler_mean"])


def test_feature_selection_is_fold_local():
    from sklearn.model_selection import StratifiedKFold
    abundance, metadata = _dataset()
    labels = pd.Categorical(metadata["Group"]).codes
    _, outer_test = next(StratifiedKFold(n_splits=5, shuffle=True, random_state=23).split(abundance.T, labels))
    abundance.loc["F_outer_test_missing"] = np.linspace(1, 2, abundance.shape[1])
    abundance.loc["F_outer_test_missing", metadata.iloc[outer_test]["Sample ID"]] = np.nan
    result = run_plsda(abundance, metadata, permutations=0, random_seed=23, max_missing_fraction=.10)
    assert "F_outer_test_missing" in result["fold_diagnostics"][0]["retained_features"]


def test_variance_filter_is_fold_local():
    from sklearn.model_selection import StratifiedKFold
    abundance, metadata = _dataset()
    labels = pd.Categorical(metadata["Group"]).codes
    _, outer_test = next(StratifiedKFold(n_splits=5, shuffle=True, random_state=24).split(abundance.T, labels))
    abundance.loc["F_outer_test_variance"] = 1.0
    abundance.loc["F_outer_test_variance", metadata.iloc[outer_test]["Sample ID"]] = 100.0
    result = run_plsda(abundance, metadata, permutations=0, random_seed=24, variance_threshold=0.0)
    assert "F_outer_test_variance" not in result["fold_diagnostics"][0]["retained_features"]


def test_preprocessing_pipeline_no_global_fit_before_cv(monkeypatch):
    abundance, metadata = _dataset()
    fitted_sample_counts = []
    original = FoldPreprocessor.fit

    def spy(self, x):
        fitted_sample_counts.append(len(x))
        return original(self, x)

    monkeypatch.setattr(FoldPreprocessor, "fit", spy)
    run_plsda(abundance, metadata, max_components=2, permutations=0, random_seed=25)
    assert fitted_sample_counts[-1] == len(metadata)
    assert all(count < len(metadata) for count in fitted_sample_counts[:-1])


def test_plsda_permutation_uses_full_nested_pipeline(monkeypatch):
    abundance, metadata = _dataset()
    original = multivariate._nested_predictions
    candidate_calls = []

    def spy(*args, **kwargs):
        candidate_calls.append(tuple(args[3]))
        return original(*args, **kwargs)

    monkeypatch.setattr(multivariate, "_nested_predictions", spy)
    result = run_plsda(abundance, metadata, max_components=3, permutations=3, random_seed=20)
    full_candidates = tuple(result["candidate_components"])
    assert candidate_calls[-3:] == [full_candidates] * 3
    assert result["permutation"]["Selected components per outer fold"].str.len().gt(0).all()


def test_plsda_oof_q2_reproducible():
    abundance, metadata = _dataset()
    first = run_plsda(abundance, metadata, permutations=2, random_seed=21)
    second = run_plsda(abundance, metadata, permutations=2, random_seed=21)
    assert first["metrics"].iloc[0]["Q2"] == second["metrics"].iloc[0]["Q2"]
    prediction_columns = [column for column in first["oof_predictions"] if column.startswith("Predicted ") and column != "Predicted class"]
    assert np.allclose(first["oof_predictions"][prediction_columns], second["oof_predictions"][prediction_columns])


def test_random_labels_not_reported_as_good_model():
    rng = np.random.default_rng(22)
    samples = [f"R{i}" for i in range(24)]
    abundance = pd.DataFrame(rng.normal(size=(20, 24)), index=[f"F{i}" for i in range(20)], columns=samples)
    labels = rng.permutation(["A"] * 12 + ["B"] * 12)
    metadata = pd.DataFrame({"Sample ID": samples, "Sample type": "biological", "Group": labels})
    result = run_plsda(abundance, metadata, permutations=2, random_seed=22)
    assert result["metrics"].iloc[0]["Predictive strength"] != "good"


def test_q2_weak_interpretation():
    interpretation = interpret_plsda_q2(0.1022)
    assert interpretation["predictive_strength"] == "weak"
    assert interpretation["label_zh"] == "预测能力较弱"


def test_significant_record_and_unique_feature_counts():
    differential = pd.DataFrame({
        "Feature ID": ["A", "B", "B", "C", "C", "C"],
        "Ion mode": ["pos"] * 6,
        "Contrast": ["C1", "C1", "C2", "C1", "C2", "C3"],
        "Significance class": ["Up", "Down", "Up", "Up", "Down", "Up"],
    })
    result = summarize_significant_features(differential)
    assert result["records"] == 6
    assert result["unique_features"] == 3
    assert result["per_contrast"] == {"C1": 3, "C2": 2, "C3": 1}


def test_manifest_fc_rule_raw():
    rule = fold_change_rule("raw", 1.2)
    assert rule["fc_formula"] == "mean(target) / mean(reference)"
    assert rule["threshold"]["lower"] == 1 / 1.2


def test_manifest_fc_rule_log2():
    rule = fold_change_rule("log2", 1.2)
    assert rule["log2fc_formula"] == "mean_log2(target) - mean_log2(reference)"
    fc, log2fc = fold_change([1], [2], "log2")
    assert fc[0] == 2 and log2fc[0] == 1


def test_manifest_fc_rule_log10():
    rule = fold_change_rule("log10", 1.2)
    assert rule["fc_formula"].startswith("10 **")
    fc, log2fc = fold_change([1], [2], "log10")
    assert fc[0] == 10 and np.isclose(log2fc[0], np.log2(10))


def test_manifest_unknown_scale_disables_fc():
    rule = fold_change_rule("unknown", 1.2)
    assert rule["computed"] is False
    assert rule["fc_formula"] is None


def test_manifest_fc_direction():
    rule = fold_change_rule("raw", 1.2)
    assert rule["contrast_direction"] == "target_over_reference"
    assert rule["contrast_label_rule"] == "TARGET vs REFERENCE"
    fc, _ = fold_change([2], [4], "raw")
    assert fc[0] == 2
