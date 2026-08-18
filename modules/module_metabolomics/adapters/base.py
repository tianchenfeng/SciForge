from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from modules.module_metabolomics.data_model import MetabolomicsDataset


OLE_MAGIC = bytes.fromhex("D0CF11E0A1B11AE1")
ZIP_MAGIC = b"PK\x03\x04"


def _looks_binary_excel(path: Path) -> bool:
    with path.open("rb") as handle:
        signature = handle.read(8)
    return signature.startswith(OLE_MAGIC) or signature.startswith(ZIP_MAGIC)


def read_table_smart(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read real Excel files and vendor TSV files with misleading .xls names."""

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input table does not exist: {path}")

    suffix = path.suffix.lower()
    if suffix in {".xlsx", ".xls"} and _looks_binary_excel(path):
        return pd.read_excel(path, **kwargs)
    if suffix == ".csv":
        return pd.read_csv(path, **kwargs)
    if suffix in {".tsv", ".txt", ".xls"}:
        return pd.read_csv(path, sep="\t", **kwargs)
    raise ValueError(f"Unsupported table format: {path.suffix}")


class MetabolomicsAdapter(ABC):
    name = "base"

    @classmethod
    @abstractmethod
    def can_load(cls, source: str | Path) -> bool:
        raise NotImplementedError

    @abstractmethod
    def load(self, source: str | Path, **kwargs) -> MetabolomicsDataset:
        raise NotImplementedError

