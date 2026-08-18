from pathlib import Path

import pandas as pd

from modules.module_metabolomics.adapters import GenericMatrixAdapter, WorkflowResultsAdapter
from modules.module_metabolomics.validation import has_blocking_errors, validate_dataset


def test_generic_matrix_adapter(tmp_path):
    matrix = pd.DataFrame({"Feature ID": ["F1", "F2"], "S1": [1, 2], "S2": [3, 4]})
    path = tmp_path / "matrix.csv"; matrix.to_csv(path, index=False)
    dataset = GenericMatrixAdapter().load(path, scale="raw")
    assert dataset.abundance.shape == (2, 2)
    assert dataset.scale == "raw"


def test_workflow_results_adapter_on_repository_fixture():
    root = Path(__file__).resolve().parents[3] / "test" / "workflow_results"
    if not root.exists():
        return
    assert WorkflowResultsAdapter.can_load(root)
    dataset = WorkflowResultsAdapter().load(root)
    assert dataset.abundance.shape == (20646, 27)
    assert len(dataset.biological_samples) == 24
    assert len(dataset.qc_samples) == 3
    assert not has_blocking_errors(validate_dataset(dataset))

