"""Optional local-data fuzzing: outputs are synthetic stress data, never study conclusions."""
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from modules.module_stats.service import run_one_way_batch_analysis, run_two_way_batch_analysis

SOURCE = Path(__file__).resolve().parents[1]/"test/凉茶渣体外发酵.xlsx"
pytestmark = pytest.mark.skipif(not SOURCE.exists(), reason="Local research fixture is intentionally not distributed")


@lru_cache
def original_numeric_values():
    values = []
    for frame in pd.read_excel(SOURCE,sheet_name=None,header=None).values():
        numeric = pd.to_numeric(pd.Series(frame.to_numpy().ravel()),errors="coerce").dropna().to_numpy(dtype=float)
        values.extend(numeric[np.isfinite(numeric)])
    assert len(values) >= 24
    return np.array(values)


@pytest.mark.parametrize("mode", ["one", "two"])
@pytest.mark.parametrize("mutation", ["baseline","large","small","outlier","missing","shuffled"])
def test_mutated_local_study_data(mode, mutation, tmp_path):
    rng = np.random.default_rng(2626)
    values = rng.choice(original_numeric_values(),size=24,replace=True)
    # Deliberately reassign groups for software robustness only.
    frame = pd.DataFrame({"Group":["C"]*12+["T"]*12,"Time":(["D1"]*6+["D2"]*6)*2,"value":values})
    if mode=="one": frame = frame.drop(columns="Time")
    if mutation=="large": frame["value"] *= 1e12
    elif mutation=="small": frame["value"] *= 1e-12
    elif mutation=="outlier": frame.loc[0,"value"] = 1e14
    elif mutation=="missing": frame.loc[[0,7,14],"value"] = np.nan
    elif mutation=="shuffled": frame = frame.sample(frac=1,random_state=26)
    source = tmp_path/"synthetic_mutated_input.csv"
    frame.to_csv(source,index=False)
    service = run_one_way_batch_analysis if mode=="one" else run_two_way_batch_analysis
    result = service(source,tmp_path/"synthetic_results")
    assert result["success"], result["message"]
    assert Path(result["excel_path"]).exists() and Path(result["word_path"]).exists()
