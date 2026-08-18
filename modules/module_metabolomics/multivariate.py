from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from modules.module_metabolomics.preprocessing import FoldPreprocessor, scale_matrix


Q2_DEFINITION = "1 - PRESS / TSS, using predictions pooled from all outer-test folds"
Q2_INTERPRETATION_NOTE = "SciForge conservative reporting rule; not a universal metabolomics industry standard."
Q2_INTERPRETATION_RULES = [
    {"minimum": None, "maximum": 0.0, "predictive_strength": "none", "label_zh": "无有效预测能力或劣于基线"},
    {"minimum": 0.0, "maximum": 0.2, "predictive_strength": "weak", "label_zh": "预测能力较弱"},
    {"minimum": 0.2, "maximum": 0.5, "predictive_strength": "limited", "label_zh": "预测能力有限，应谨慎解释"},
    {"minimum": 0.5, "maximum": None, "predictive_strength": "good", "label_zh": "在当前交叉验证设置下具有较好的预测表现，但仍需结合置换检验和样本量判断"},
]


def interpret_plsda_q2(q2: float) -> dict:
    value = float(q2)
    if not np.isfinite(value):
        return {"predictive_strength": "unavailable", "label_zh": "预测能力不可用", "q2": None}
    if value < 0:
        rule = Q2_INTERPRETATION_RULES[0]
    elif value < 0.2:
        rule = Q2_INTERPRETATION_RULES[1]
    elif value < 0.5:
        rule = Q2_INTERPRETATION_RULES[2]
    else:
        rule = Q2_INTERPRETATION_RULES[3]
    return {"predictive_strength": rule["predictive_strength"], "label_zh": rule["label_zh"], "q2": value}


def _rows(x, indices):
    return x.iloc[indices] if isinstance(x, pd.DataFrame) else x[indices]


def run_pca(abundance: pd.DataFrame, metadata: pd.DataFrame, scaling: str = "pareto", n_components: int = 10) -> dict:
    from sklearn.decomposition import PCA
    sample_ids = [sample for sample in metadata["Sample ID"].astype(str) if sample in abundance.columns]
    x = abundance[sample_ids].T.apply(pd.to_numeric, errors="coerce")
    x = x.fillna(x.median(axis=0)).fillna(0)
    scaled = scale_matrix(x, scaling)
    components = min(n_components, scaled.shape[0] - 1, scaled.shape[1])
    if components < 2:
        raise ValueError("PCA requires at least three samples and two variable dimensions.")
    model = PCA(n_components=components, svd_solver="full")
    coordinates = model.fit_transform(scaled)
    names = [f"PC{i + 1}" for i in range(components)]
    scores = pd.DataFrame(coordinates, columns=names)
    scores.insert(0, "Sample ID", sample_ids)
    scores = scores.merge(metadata, on="Sample ID", how="left")
    loadings = pd.DataFrame(model.components_.T, index=scaled.columns, columns=names).reset_index(names="Feature ID")
    variance = pd.DataFrame({"Component": names, "Explained variance": model.explained_variance_ratio_, "Cumulative variance": np.cumsum(model.explained_variance_ratio_)})
    return {"scores": scores, "loadings": loadings, "variance": variance, "scaling": scaling, "model": model}


def _fit_pls(x, y, components):
    from sklearn.cross_decomposition import PLSRegression
    model = PLSRegression(n_components=components, scale=False, max_iter=1000)
    model.fit(x, y)
    return model


def _q2(y_true, y_pred) -> float:
    press = float(np.sum((y_true - y_pred) ** 2))
    centered = float(np.sum((y_true - y_true.mean(axis=0)) ** 2))
    return 1 - press / centered if centered > 0 else np.nan


def _balanced_accuracy(y_true, prediction) -> float:
    from sklearn.metrics import balanced_accuracy_score
    return float(balanced_accuracy_score(np.argmax(y_true, axis=1), np.argmax(prediction, axis=1)))


def _inner_select(x, y, labels, candidates, seed, scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold, return_diagnostics=False):
    from sklearn.model_selection import StratifiedKFold
    counts = np.bincount(labels)
    splits = min(3, int(counts.min()))
    if splits < 2:
        return (1, []) if return_diagnostics else 1
    cv = StratifiedKFold(n_splits=splits, shuffle=True, random_state=seed)
    best = (float("-inf"), 1)
    diagnostics = []
    predictions = {components: np.zeros_like(y, dtype=float) for components in candidates}
    valid = {components: True for components in candidates}
    for train, test in cv.split(x, labels):
        processor = FoldPreprocessor(scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold)
        x_train = processor.fit_transform(_rows(x, train))
        x_test = processor.transform(_rows(x, test))
        max_allowed = min(len(train) - 1, x_train.shape[1], y.shape[1] + 2)
        for components in candidates:
            if components > max_allowed:
                valid[components] = False
                continue
            predictions[components][test] = _fit_pls(x_train, y[train], components).predict(x_test)
        if return_diagnostics:
            diagnostics.append({"components_evaluated": [component for component in candidates if component <= max_allowed], "train_indices": train.copy(), "validation_indices": test.copy(), **processor.diagnostics()})
    for components in candidates:
        if valid[components]:
            score = _q2(y, predictions[components])
            if score > best[0]: best = (score, components)
    return (best[1], diagnostics) if return_diagnostics else best[1]


def _nested_predictions(x, y, labels, candidates, folds, seed, scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold, collect_diagnostics=False):
    from sklearn.model_selection import StratifiedKFold
    min_group = int(np.bincount(labels).min())
    splits = min(folds, min_group)
    if splits < 2:
        raise ValueError("PLS-DA cross-validation requires at least two samples in every group.")
    cv = StratifiedKFold(n_splits=splits, shuffle=True, random_state=seed)
    prediction = np.zeros_like(y, dtype=float); selected = []; diagnostics = []; fold_ids = np.zeros(len(labels), dtype=int)
    for fold, (train, test) in enumerate(cv.split(x, labels)):
        outer_train = _rows(x, train)
        inner = _inner_select(outer_train, y[train], labels[train], candidates, seed + fold + 1, scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold, collect_diagnostics)
        components, inner_diagnostics = inner if collect_diagnostics else (inner, [])
        selected.append(components)
        processor = FoldPreprocessor(scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold)
        x_train = processor.fit_transform(outer_train); x_test = processor.transform(_rows(x, test))
        prediction[test] = _fit_pls(x_train, y[train], components).predict(x_test)
        fold_ids[test] = fold + 1
        if collect_diagnostics:
            diagnostics.append({"outer_fold": fold + 1, "train_indices": train.copy(), "test_indices": test.copy(), "selected_components": components, "inner": inner_diagnostics, **processor.diagnostics()})
    return prediction, selected, diagnostics, fold_ids


def _vip(model, feature_names):
    t = model.x_scores_; w = model.x_weights_; q = model.y_loadings_
    p, h = w.shape
    s = np.diag(t.T @ t @ q.T @ q)
    total = float(s.sum())
    if total <= 0:
        values = np.full(p, np.nan)
    else:
        norms = np.linalg.norm(w, axis=0)
        normalized = w / np.where(norms == 0, 1, norms)
        values = np.sqrt(p * ((normalized ** 2) @ s) / total)
    return pd.DataFrame({"Feature ID": feature_names, "VIP": values}).sort_values("VIP", ascending=False)


def run_plsda(
    abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    scaling: str = "pareto",
    max_components: int = 5,
    cv_folds: int = 5,
    permutations: int = 200,
    random_seed: int = 2026,
    progress=None,
    missing_method: str = "none",
    transformation: str = "none",
    input_scale: str = "raw",
    max_missing_fraction: float = 1.0,
    variance_threshold: float | None = None,
) -> dict:
    meta = metadata[metadata["Sample type"].astype(str).str.lower().ne("qc")].copy()
    sample_ids = [sample for sample in meta["Sample ID"].astype(str) if sample in abundance.columns]
    meta = meta.set_index("Sample ID").loc[sample_ids].reset_index()
    groups = meta["Group"].astype(str)
    if groups.nunique() < 2:
        raise ValueError("PLS-DA requires at least two groups.")
    x_df = abundance[sample_ids].T.apply(pd.to_numeric, errors="coerce")
    classes = list(pd.unique(groups)); labels = pd.Categorical(groups, categories=classes).codes
    y = pd.get_dummies(groups).reindex(columns=classes, fill_value=0).to_numpy(dtype=float)
    x = x_df
    upper = min(max_components, x.shape[0] - 2, x.shape[1], max(1, y.shape[1] + 2))
    candidates = list(range(1, max(upper, 1) + 1))
    prediction, selected, diagnostics, fold_ids = _nested_predictions(x, y, labels, candidates, cv_folds, random_seed, scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold, True)
    selected_components = Counter(selected).most_common(1)[0][0]
    q2 = _q2(y, prediction); accuracy = _balanced_accuracy(y, prediction)

    cv_rows = []
    for components in candidates:
        pred, _, _, _ = _nested_predictions(x, y, labels, [components], cv_folds, random_seed, scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold)
        cv_rows.append({"Components": components, "Q2": _q2(y, pred), "Balanced accuracy": _balanced_accuracy(y, pred)})
    cv_table = pd.DataFrame(cv_rows)
    full_processor = FoldPreprocessor(scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold)
    x_full = full_processor.fit_transform(x)
    model = _fit_pls(x_full, y, selected_components)
    fitted = model.predict(x_full)
    r2y = 1 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean(axis=0)) ** 2)
    reconstructed = model.x_scores_ @ model.x_loadings_.T
    r2x = 1 - np.sum((x_full - reconstructed) ** 2) / np.sum((x_full - x_full.mean(axis=0)) ** 2)
    score_columns = [f"PLS{i + 1}" for i in range(model.x_scores_.shape[1])]
    scores = pd.DataFrame(model.x_scores_, columns=score_columns)
    scores.insert(0, "Sample ID", sample_ids); scores = scores.merge(meta, on="Sample ID", how="left")
    vip = _vip(model, full_processor.retained_features_)
    oof = pd.DataFrame(prediction, columns=[f"Predicted {name}" for name in classes])
    oof.insert(0, "Outer fold", fold_ids); oof.insert(0, "True class", groups.to_numpy()); oof.insert(0, "Sample ID", sample_ids)
    oof["Predicted class"] = [classes[index] for index in np.argmax(prediction, axis=1)]

    rng = np.random.default_rng(random_seed); permutation_rows = []
    for index in range(permutations):
        permuted_labels = rng.permutation(labels)
        permuted_y = np.eye(len(classes))[permuted_labels]
        pred, perm_selected, _, _ = _nested_predictions(x, permuted_y, permuted_labels, candidates, cv_folds, random_seed + index + 1000, scaling, missing_method, transformation, input_scale, max_missing_fraction, variance_threshold)
        permutation_rows.append({"Permutation": index + 1, "Q2": _q2(permuted_y, pred), "Balanced accuracy": _balanced_accuracy(permuted_y, pred), "Selected components per outer fold": ",".join(map(str, perm_selected))})
        if callable(progress) and (index == permutations - 1 or (index + 1) % max(permutations // 100, 1) == 0): progress((index + 1) / permutations)
    permutation_table = pd.DataFrame(permutation_rows)
    permutation_p = (1 + int((permutation_table["Q2"] >= q2).sum())) / (permutations + 1) if permutations else np.nan
    interpretation = interpret_plsda_q2(q2)
    outer_folds = min(cv_folds, int(np.bincount(labels).min()))
    inner_folds = min(3, int(np.bincount(labels).min()) - 1)
    metrics = pd.DataFrame([{"Components": selected_components, "R2X": r2x, "R2Y": r2y, "Q2": q2, "Q2 definition": Q2_DEFINITION, "Predictive strength": interpretation["predictive_strength"], "Predictive interpretation": interpretation["label_zh"], "Balanced accuracy": accuracy, "Outer CV folds": outer_folds, "Outer repeats": 1, "Inner CV folds": inner_folds, "Candidate components": ",".join(map(str, candidates)), "Selected components per outer fold": ",".join(map(str, selected)), "Retained features per outer fold": ",".join(str(len(item["retained_features"])) for item in diagnostics), "Scaling": scaling, "Feature selection fitted within folds": True, "Scaling fitted within folds": True, "Maximum missing fraction": max_missing_fraction, "Variance threshold": variance_threshold, "Permutations": permutations, "Permutation P": permutation_p, "Permutation status": "passed" if np.isfinite(permutation_p) and permutation_p < 0.05 and q2 > 0 else "not passed", "Random seed": random_seed}])
    return {"scores": scores, "vip": vip, "cv": cv_table, "oof_predictions": oof, "permutation": permutation_table, "metrics": metrics, "model": model, "full_preprocessor": full_processor, "fold_diagnostics": diagnostics, "scaling": scaling, "classes": classes, "candidate_components": candidates, "selected_components_per_outer_fold": selected}
