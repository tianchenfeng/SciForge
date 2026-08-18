from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


SAMPLE_COLUMNS = ["Sample ID", "Sample type", "Group", "Factor A", "Factor B"]


@dataclass
class MetabolomicsDataset:
    """Vendor-neutral, feature-by-sample metabolomics dataset."""

    abundance: pd.DataFrame
    samples: pd.DataFrame
    features: pd.DataFrame
    modeling_abundance: pd.DataFrame | None = None
    modeling_scale: str | None = None
    raw_abundance: pd.DataFrame | None = None
    qc_rsd: pd.DataFrame | None = None
    provenance: dict[str, Any] = field(default_factory=dict)
    vendor_results: dict[str, Any] = field(default_factory=dict)
    scale: str = "unknown"
    normalized: bool | None = None
    source: Path | None = None

    def copy(self) -> "MetabolomicsDataset":
        return MetabolomicsDataset(
            abundance=self.abundance.copy(),
            samples=self.samples.copy(),
            features=self.features.copy(),
            modeling_abundance=None if self.modeling_abundance is None else self.modeling_abundance.copy(),
            modeling_scale=self.modeling_scale,
            raw_abundance=None if self.raw_abundance is None else self.raw_abundance.copy(),
            qc_rsd=None if self.qc_rsd is None else self.qc_rsd.copy(),
            provenance=dict(self.provenance),
            vendor_results=dict(self.vendor_results),
            scale=self.scale,
            normalized=self.normalized,
            source=self.source,
        )

    @property
    def biological_samples(self) -> list[str]:
        mask = self.samples["Sample type"].astype(str).str.lower().ne("qc")
        return self.samples.loc[mask, "Sample ID"].astype(str).tolist()

    @property
    def qc_samples(self) -> list[str]:
        mask = self.samples["Sample type"].astype(str).str.lower().eq("qc")
        return self.samples.loc[mask, "Sample ID"].astype(str).tolist()

    def freeze_metadata(self) -> list[dict[str, Any]]:
        return self.samples[SAMPLE_COLUMNS].fillna("").to_dict("records")
