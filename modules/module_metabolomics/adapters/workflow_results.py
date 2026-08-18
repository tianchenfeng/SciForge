from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from modules.module_metabolomics.adapters.base import MetabolomicsAdapter, read_table_smart
from modules.module_metabolomics.data_model import MetabolomicsDataset


FEATURE_COLUMN_MAP = {
    "metab_id": "Feature ID",
    "Metabolite": "Metabolite name",
    "metab": "Metabolite name",
    "Mode": "Ion mode",
    "mode": "Ion mode",
    "m/z": "m/z",
    "Retention time": "RT",
    "Formula": "Formula",
    "formula": "Formula",
    "KEGG Compound ID": "KEGG",
    "compound_id": "KEGG",
    "pathway_id": "Pathway IDs",
    "description": "Pathway names",
    "hmdb_id": "HMDB",
    "Library ID": "Library ID",
    "level": "Annotation level",
    "Fragmentation Score": "MS/MS score",
    "Adducts": "Adduct",
}


class WorkflowResultsAdapter(MetabolomicsAdapter):
    name = "workflow_results"

    @classmethod
    def can_load(cls, source: str | Path) -> bool:
        root = Path(source)
        return (root / "1.Preprocess" / "mix" / "metab_abund.txt").exists()

    @staticmethod
    def _infer_sample(sample_id: str) -> dict[str, str]:
        sample_id = str(sample_id).strip()
        if re.match(r"^QC", sample_id, flags=re.I):
            return {"Sample ID": sample_id, "Sample type": "QC", "Group": "QC", "Factor A": "", "Factor B": ""}
        match = re.match(r"^(CON|LP|SC|MIX)", sample_id, flags=re.I)
        group = match.group(1).upper() if match else ""
        return {"Sample ID": sample_id, "Sample type": "biological", "Group": group, "Factor A": group, "Factor B": ""}

    @staticmethod
    def _load_abundance(path: Path) -> pd.DataFrame:
        table = read_table_smart(path, low_memory=False)
        ids = table.iloc[:, 0].astype(str).str.strip()
        values = table.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        values.index = ids
        values.index.name = "Feature ID"
        values.columns = values.columns.astype(str).str.strip()
        return values

    @staticmethod
    def _load_rsd(path: Path) -> pd.DataFrame | None:
        if not path.exists():
            return None
        table = pd.read_csv(path, sep="\t", header=None, names=["Feature ID", "QC RSD"], usecols=[0, 1])
        if not table.empty and str(table.iloc[0, 0]).strip() in {"0", "Feature ID", "metab_id"}:
            table = table.iloc[1:].copy()
        table["Feature ID"] = table["Feature ID"].astype(str).str.strip()
        table["QC RSD"] = pd.to_numeric(table["QC RSD"], errors="coerce")
        return table.dropna(subset=["Feature ID"]).drop_duplicates("Feature ID", keep="first")

    @staticmethod
    def _load_features(root: Path) -> pd.DataFrame:
        desc = read_table_smart(root / "1.Preprocess" / "mix" / "metab_desc.txt", low_memory=False)
        keep = [c for c in FEATURE_COLUMN_MAP if c in desc.columns]
        features = desc[keep].rename(columns=FEATURE_COLUMN_MAP).copy()
        features["Feature ID"] = features["Feature ID"].astype(str).str.strip()

        anno_path = root / "3.Anno" / "04.AnnoOverview" / "Anno_data.txt"
        if anno_path.exists():
            anno = read_table_smart(anno_path, low_memory=False)
            keep_anno = [c for c in FEATURE_COLUMN_MAP if c in anno.columns]
            anno = anno[keep_anno].rename(columns=FEATURE_COLUMN_MAP)
            anno["Feature ID"] = anno["Feature ID"].astype(str).str.strip()
            anno = anno.drop_duplicates("Feature ID")
            features = features.set_index("Feature ID")
            anno = anno.set_index("Feature ID")
            for column in anno.columns:
                incoming = anno[column].replace("-", pd.NA)
                if column not in features:
                    features[column] = incoming
                else:
                    current = features[column].replace(["-", "unknown", ""], pd.NA)
                    features[column] = incoming.combine_first(current)
            features = features.reset_index()

        for column in ["Ion mode", "Metabolite name", "m/z", "RT", "Formula", "KEGG", "HMDB", "Library ID", "Annotation level", "MS/MS score", "Pathway IDs", "Pathway names"]:
            if column not in features:
                features[column] = ""
        features["Identified"] = ~features["Annotation level"].fillna("").astype(str).isin(["", "-", "unknown"])
        features.loc[~features["Identified"], "Metabolite name"] = ""
        return features

    @staticmethod
    def _vendor_result_paths(root: Path) -> dict[str, str]:
        patterns = {
            "pca_model": "2.SampleComp/02.ExpPCA/mix/PCA.model.xls",
            "pca_scores": "2.SampleComp/02.ExpPCA/mix/PCA.sites.xls",
            "plsda_model": "2.SampleComp/04.ExpPLSDA/mix/PLS-DA.model.xls",
            "plsda_permutation": "2.SampleComp/04.ExpPLSDA/mix/PLS-DA.permMN.xls",
            "oplsda_root": "4.ExpDiff/01.TwoGroupExpDiff/mix/DiffMulStat",
            "vendor_differential_root": "4.ExpDiff/01.TwoGroupExpDiff/mix/DiffTest",
            "vendor_pathway_root": "5.Metabset/07.MetabsetEnrich",
        }
        return {key: str(root / value) for key, value in patterns.items() if (root / value).exists()}

    def load(self, source: str | Path, **kwargs) -> MetabolomicsDataset:
        root = Path(source)
        if not self.can_load(root):
            raise ValueError(f"Not a supported workflow_results directory: {root}")
        abundance = self._load_abundance(root / "1.Preprocess" / "mix" / "metab_abund.txt")
        raw_path = root / "1.Preprocess" / "org_mix" / "metab_abund.txt"
        raw = self._load_abundance(raw_path) if raw_path.exists() else None
        samples = pd.DataFrame([self._infer_sample(column) for column in abundance.columns])
        features = self._load_features(root)
        rsd = self._load_rsd(root / "1.Preprocess" / "mix" / "metab_abund.txt_rsd.txt")
        return MetabolomicsDataset(
            abundance=abundance,
            raw_abundance=raw,
            samples=samples,
            features=features,
            qc_rsd=rsd,
            provenance={
                "adapter": self.name,
                "input": str(root),
                "abundance": str(root / "1.Preprocess" / "mix" / "metab_abund.txt"),
                "raw_abundance": str(raw_path) if raw_path.exists() else None,
                "scale_note": "Vendor processed scale could not be proven from the supplied files.",
            },
            vendor_results=self._vendor_result_paths(root),
            scale="unknown",
            normalized=True,
            source=root,
        )
