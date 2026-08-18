from __future__ import annotations

from pathlib import Path

import pandas as pd

from modules.module_metabolomics.adapters.base import MetabolomicsAdapter, read_table_smart
from modules.module_metabolomics.data_model import MetabolomicsDataset, SAMPLE_COLUMNS


class GenericMatrixAdapter(MetabolomicsAdapter):
    name = "generic_matrix"

    @classmethod
    def can_load(cls, source: str | Path) -> bool:
        path = Path(source)
        return path.is_file() and path.suffix.lower() in {".xlsx", ".xls", ".csv", ".tsv", ".txt"}

    def load(
        self,
        source: str | Path,
        metadata_file: str | Path | None = None,
        annotation_file: str | Path | None = None,
        scale: str = "unknown",
        normalized: bool | None = None,
    ) -> MetabolomicsDataset:
        source = Path(source)
        table = read_table_smart(source, low_memory=False)
        if table.shape[1] < 3:
            raise ValueError("Generic abundance matrix requires a feature ID column and at least two samples.")
        feature_id_col = str(table.columns[0])
        feature_ids = table.iloc[:, 0].astype(str).str.strip()
        abundance = table.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        abundance.index = feature_ids
        abundance.index.name = "Feature ID"
        abundance.columns = abundance.columns.astype(str).str.strip()

        if metadata_file:
            samples = read_table_smart(metadata_file, low_memory=False)
            lower = {str(c).strip().lower(): c for c in samples.columns}
            sample_col = lower.get("sample id") or lower.get("sample_id") or lower.get("sample")
            if sample_col is None:
                raise ValueError("Sample metadata requires a Sample ID column.")
            samples = samples.rename(columns={sample_col: "Sample ID"})
        else:
            samples = pd.DataFrame({"Sample ID": abundance.columns})
            samples["Sample type"] = "biological"
            samples["Group"] = ""
            samples["Factor A"] = ""
            samples["Factor B"] = ""

        for col in SAMPLE_COLUMNS:
            if col not in samples:
                samples[col] = ""
        samples = samples[SAMPLE_COLUMNS].copy()
        samples["Sample ID"] = samples["Sample ID"].astype(str).str.strip()

        features = pd.DataFrame({"Feature ID": feature_ids})
        features["Ion mode"] = "unknown"
        features["Metabolite name"] = ""
        features["Annotation level"] = ""
        if annotation_file:
            anno = read_table_smart(annotation_file, low_memory=False)
            anno = anno.rename(columns={anno.columns[0]: "Feature ID"})
            anno["Feature ID"] = anno["Feature ID"].astype(str).str.strip()
            features = features.merge(anno, on="Feature ID", how="left", suffixes=("", "_input"))

        return MetabolomicsDataset(
            abundance=abundance,
            samples=samples,
            features=features,
            provenance={"adapter": self.name, "input": str(source), "feature_id_column": feature_id_col},
            scale=scale,
            normalized=normalized,
            source=source,
        )

