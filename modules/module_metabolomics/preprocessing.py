from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

from modules.module_metabolomics.data_model import MetabolomicsDataset


@dataclass
class PreprocessingConfig:
    missing_method: str = "none"
    transformation: str = "none"
    scaling: str = "pareto"
    normalization: str = "existing"
    max_missing_fraction: float = 1.0
    qc_rsd_threshold: float | None = None
    confirmed_input_scale: str = "unknown"

    def to_dict(self) -> dict:
        return asdict(self)


class MetabolomicsScaler:
    """Feature-wise scaler whose parameters are fitted on training samples only."""

    def __init__(self, method: str = "pareto"):
        self.method = str(method).lower().replace(" ", "_")
        self.mean_: np.ndarray | None = None
        self.std_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None
        self.n_samples_seen_: int = 0

    def fit(self, x):
        values = np.asarray(x, dtype=float)
        if values.ndim != 2:
            raise ValueError("MetabolomicsScaler expects a two-dimensional matrix.")
        self.n_samples_seen_ = int(values.shape[0])
        self.mean_ = np.mean(values, axis=0)
        std = np.std(values, axis=0, ddof=1)
        std = np.where(np.isfinite(std) & (std > 0), std, 1.0)
        self.std_ = std
        if self.method in {"none", ""}:
            self.scale_ = np.ones(values.shape[1], dtype=float)
        elif self.method in {"mean_centering", "center", "mean_center"}:
            self.scale_ = np.ones(values.shape[1], dtype=float)
        elif self.method in {"autoscaling", "uv", "unit_variance"}:
            self.scale_ = std
        elif self.method == "pareto":
            self.scale_ = np.sqrt(std)
        else:
            raise ValueError(f"Unknown scaling method: {self.method}")
        return self

    def transform(self, x):
        if self.mean_ is None or self.scale_ is None:
            raise RuntimeError("MetabolomicsScaler must be fitted before transform.")
        values = np.asarray(x, dtype=float)
        if self.method in {"none", ""}:
            return values.copy()
        return (values - self.mean_) / self.scale_

    def fit_transform(self, x):
        return self.fit(x).transform(x)


class FoldPreprocessor:
    """Training-fold-only feature selection, imputation, transformation, and scaling."""

    def __init__(
        self,
        scaling: str = "pareto",
        missing_method: str = "none",
        transformation: str = "none",
        input_scale: str = "raw",
        max_missing_fraction: float = 1.0,
        variance_threshold: float | None = None,
    ):
        self.scaling = scaling
        self.missing_method = str(missing_method).lower().replace("-", "_").replace(" ", "_")
        self.transformation = str(transformation).lower()
        self.input_scale = str(input_scale).lower()
        self.max_missing_fraction = float(max_missing_fraction)
        self.variance_threshold = variance_threshold
        self.feature_names_in_: list[str] | None = None
        self.missing_fraction_: np.ndarray | None = None
        self.missing_threshold_: float = self.max_missing_fraction
        self.missing_retained_features_: list[str] | None = None
        self.retained_features_: list[str] | None = None
        self.imputation_values_: np.ndarray | None = None
        self.variance_: np.ndarray | None = None
        self.pseudocount_: float | None = None
        self.knn_ = None
        self.scaler_ = MetabolomicsScaler(scaling)

    def _frame(self, x, fitting=False) -> pd.DataFrame:
        if isinstance(x, pd.DataFrame):
            frame = x if all(pd.api.types.is_numeric_dtype(dtype) for dtype in x.dtypes) else x.apply(pd.to_numeric, errors="coerce")
            if fitting:
                self.feature_names_in_ = list(map(str, frame.columns))
                frame.columns = self.feature_names_in_
                return frame
            return frame.reindex(columns=self.feature_names_in_)
        values = np.asarray(x, dtype=float)
        if values.ndim != 2:
            raise ValueError("FoldPreprocessor expects a two-dimensional matrix.")
        if fitting:
            self.feature_names_in_ = [str(index) for index in range(values.shape[1])]
        if self.feature_names_in_ is None or values.shape[1] != len(self.feature_names_in_):
            raise ValueError("FoldPreprocessor transform matrix does not match fitted features.")
        return pd.DataFrame(values, columns=self.feature_names_in_)

    def _fit_imputer(self, values: np.ndarray) -> np.ndarray:
        if self.missing_method == "knn":
            from sklearn.impute import KNNImputer
            self.knn_ = KNNImputer(n_neighbors=min(5, max(len(values) - 1, 1)), keep_empty_features=True).fit(values)
            self.imputation_values_ = np.asarray(self.knn_._fit_X, dtype=float).copy()
            return self.knn_.transform(values)
        if self.missing_method in {"half_minimum", "half_min", "group_wise_half_minimum", "group_half_minimum", "group_half_min"}:
            positive = np.where(values > 0, values, np.nan)
            with np.errstate(all="ignore"):
                fill = np.nanmin(positive, axis=0) / 2.0
                global_fill = np.nanmin(positive) / 2.0
            if not np.isfinite(global_fill) or global_fill <= 0:
                global_fill = 0.0
            self.imputation_values_ = np.where(np.isfinite(fill), fill, global_fill)
        else:
            with np.errstate(all="ignore"):
                fill = np.nanmedian(values, axis=0)
            self.imputation_values_ = np.where(np.isfinite(fill), fill, 0.0)
        return np.where(np.isfinite(values), values, self.imputation_values_)

    def _impute(self, values: np.ndarray) -> np.ndarray:
        if self.knn_ is not None:
            return self.knn_.transform(values)
        return np.where(np.isfinite(values), values, self.imputation_values_)

    def _fit_transform_values(self, values: np.ndarray) -> np.ndarray:
        imputed = self._fit_imputer(values)
        if self.transformation in {"log2", "log10"}:
            if self.input_scale in {"log2", "log10"}:
                raise ValueError(f"Input is already declared as {self.input_scale}; another log transform is not allowed.")
            positive = imputed[imputed > 0]
            self.pseudocount_ = float(positive.min() / 2.0) if positive.size else None
            if self.pseudocount_ is None or not np.isfinite(self.pseudocount_) or self.pseudocount_ <= 0:
                raise ValueError("Log transformation requires a positive value in every training fold.")
        return self._transform_values(imputed)

    def _transform_values(self, values: np.ndarray) -> np.ndarray:
        if self.transformation not in {"log2", "log10"}:
            return values
        safe = np.where(values > 0, values, self.pseudocount_)
        return np.log2(safe) if self.transformation == "log2" else np.log10(safe)

    def fit(self, x):
        frame = self._frame(x, fitting=True)
        self.missing_fraction_ = frame.isna().mean(axis=0).to_numpy(dtype=float)
        missing_keep = self.missing_fraction_ <= self.missing_threshold_
        if not missing_keep.any():
            raise ValueError("No features remain after training-fold missing-fraction filtering.")
        self.missing_retained_features_ = list(frame.columns[missing_keep])
        transformed = self._fit_transform_values(frame.loc[:, self.missing_retained_features_].to_numpy(dtype=float))
        self.variance_ = np.var(transformed, axis=0, ddof=1)
        variance_keep = np.ones(transformed.shape[1], dtype=bool)
        if self.variance_threshold is not None:
            variance_keep = np.isfinite(self.variance_) & (self.variance_ > float(self.variance_threshold))
        if not variance_keep.any():
            raise ValueError("No features remain after training-fold variance filtering.")
        self._variance_keep_ = variance_keep
        self.retained_features_ = [name for name, keep in zip(self.missing_retained_features_, variance_keep) if keep]
        self.scaler_.fit(transformed[:, variance_keep])
        return self

    def transform(self, x) -> np.ndarray:
        if self.retained_features_ is None or self.missing_retained_features_ is None:
            raise RuntimeError("FoldPreprocessor must be fitted before transform.")
        frame = self._frame(x)
        values = frame.loc[:, self.missing_retained_features_].to_numpy(dtype=float)
        transformed = self._transform_values(self._impute(values))[:, self._variance_keep_]
        return self.scaler_.transform(transformed)

    def fit_transform(self, x) -> np.ndarray:
        return self.fit(x).transform(x)

    def diagnostics(self) -> dict:
        return {
            "n_samples_fitted": self.scaler_.n_samples_seen_,
            "retained_features": list(self.retained_features_ or []),
            "missing_threshold": self.missing_threshold_,
            "missing_fraction": None if self.missing_fraction_ is None else self.missing_fraction_.copy(),
            "imputation_values": None if self.imputation_values_ is None else self.imputation_values_.copy(),
            "variance_threshold": self.variance_threshold,
            "variance": None if self.variance_ is None else self.variance_.copy(),
            "scaler_mean": self.scaler_.mean_.copy(),
            "scaler_std": self.scaler_.std_.copy(),
            "pseudocount": self.pseudocount_,
        }


def _half_minimum(values: pd.DataFrame) -> pd.DataFrame:
    result = values.copy()
    positive = result.where(result > 0)
    fill = positive.min(axis=1) / 2.0
    global_fill = positive.min().min() / 2.0
    fill = fill.fillna(global_fill)
    return result.T.fillna(fill).T


def _group_half_minimum(values: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    result = values.copy()
    global_fill = result.where(result > 0).min().min() / 2.0
    meta = samples.set_index("Sample ID")
    for group, ids in meta[meta["Sample type"].str.lower().ne("qc")].groupby("Group").groups.items():
        columns = [sample for sample in ids if sample in result.columns]
        if not columns:
            continue
        block = result[columns]
        fill = (block.where(block > 0).min(axis=1) / 2.0).fillna(global_fill)
        result.loc[:, columns] = block.T.fillna(fill).T
    qc_columns = [c for c in result if c in set(meta[meta["Sample type"].str.lower().eq("qc")].index)]
    if qc_columns:
        block = result[qc_columns]
        fill = (block.where(block > 0).min(axis=1) / 2.0).fillna(global_fill)
        result.loc[:, qc_columns] = block.T.fillna(fill).T
    return result


def impute_missing(values: pd.DataFrame, method: str, samples: pd.DataFrame) -> pd.DataFrame:
    method = method.lower().replace("-", "_").replace(" ", "_")
    if method in {"none", "no", ""}:
        return values.copy()
    if method in {"half_minimum", "half_min"}:
        return _half_minimum(values)
    if method in {"group_wise_half_minimum", "group_half_minimum", "group_half_min"}:
        return _group_half_minimum(values, samples)
    if method == "median":
        return values.T.fillna(values.median(axis=1)).T
    if method == "knn":
        try:
            from sklearn.impute import KNNImputer
        except ImportError as exc:
            raise RuntimeError("KNN imputation requires scikit-learn.") from exc
        imputer = KNNImputer(n_neighbors=min(5, max(values.shape[1] - 1, 1)))
        return pd.DataFrame(imputer.fit_transform(values.T).T, index=values.index, columns=values.columns)
    raise ValueError(f"Unknown missing-value method: {method}")


def transform_abundance(values: pd.DataFrame, transformation: str, current_scale: str) -> tuple[pd.DataFrame, str]:
    transformation = transformation.lower()
    if transformation in {"none", ""}:
        return values.copy(), current_scale
    if current_scale in {"log2", "log10"}:
        raise ValueError(f"Input is already declared as {current_scale}; another log transform is not allowed.")
    if (values < 0).to_numpy().any():
        raise ValueError("Log transformation cannot be applied to negative values.")
    positive = values.where(values > 0)
    half_global_min = positive.min().min() / 2.0
    if not np.isfinite(half_global_min) or half_global_min <= 0:
        raise ValueError("Log transformation requires at least one positive value.")
    safe = values.mask(values <= 0, half_global_min)
    if transformation == "log2":
        return np.log2(safe), "log2"
    if transformation == "log10":
        return np.log10(safe), "log10"
    raise ValueError(f"Unknown transformation: {transformation}")


def scale_matrix(sample_by_feature: pd.DataFrame, method: str) -> pd.DataFrame:
    scaler = MetabolomicsScaler(method)
    scaled = scaler.fit_transform(sample_by_feature.to_numpy(dtype=float))
    return pd.DataFrame(scaled, index=sample_by_feature.index, columns=sample_by_feature.columns)


def preprocess_dataset(dataset: MetabolomicsDataset, config: PreprocessingConfig) -> tuple[MetabolomicsDataset, pd.DataFrame]:
    result = dataset.copy()
    values = result.abundance.apply(pd.to_numeric, errors="coerce")
    rows = [{"Stage": "Raw", "Features": len(values)}]

    declared_scale = config.confirmed_input_scale if config.confirmed_input_scale != "unknown" else result.scale
    # PLS-DA consumes this pre-selection matrix and fits all data-dependent
    # feature selection and preprocessing independently inside each CV fold.
    result.modeling_abundance = values.copy()
    result.modeling_scale = declared_scale

    keep_missing = values.isna().mean(axis=1) <= config.max_missing_fraction
    values = values.loc[keep_missing]
    rows.append({"Stage": "Missing-filtered", "Features": len(values)})

    if config.qc_rsd_threshold is not None and result.qc_rsd is not None:
        rsd = result.qc_rsd.set_index("Feature ID")["QC RSD"]
        keep_ids = rsd[rsd <= config.qc_rsd_threshold].index
        values = values.loc[values.index.intersection(keep_ids)]
    rows.append({"Stage": "QC-filtered", "Features": len(values)})

    values = impute_missing(values, config.missing_method, result.samples)
    transformation_pseudocount = None
    if config.transformation.lower() in {"log2", "log10"}:
        positive = values.where(values > 0)
        candidate = positive.min().min() / 2.0
        transformation_pseudocount = float(candidate) if np.isfinite(candidate) and candidate > 0 else None
    values, output_scale = transform_abundance(values, config.transformation, declared_scale)
    result.abundance = values
    result.features = result.features[result.features["Feature ID"].astype(str).isin(values.index.astype(str))].copy()
    result.scale = output_scale
    result.provenance["preprocessing"] = {**config.to_dict(), "transformation_pseudocount": transformation_pseudocount}
    rows.append({"Stage": "Final", "Features": len(values)})
    return result, pd.DataFrame(rows)
