"""Run every bundled 16S stage with a small synthetic study, in both supported modes."""
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from modules.module_16s.config import PIPELINE_ORDER
from modules.module_16s.service import SixteenSService

ROOT = Path(__file__).resolve().parents[1]


def make_input(parent):
    rng = np.random.default_rng(2626)
    source = parent / "input"
    qiime = source / "qiime2"
    picrust = source / "picrust2"
    qiime.mkdir(parents=True)
    picrust.mkdir()
    ids = [f"S{i}" for i in range(12)]
    pd.DataFrame({"SampleID":ids,"Treatment":["C"]*6+["T"]*6,
                  "Time":[3]*3+[7]*3+[3]*3+[7]*3}).to_csv(qiime/"metadata.tsv",sep="\t",index=False)
    features = [f"F{i}" for i in range(32)]
    counts = pd.DataFrame(rng.integers(50,500,size=(32,12)),index=features,columns=ids)
    counts.index.name = "#OTU ID"
    (qiime/"feature-table.tsv").write_text("# Constructed from biom file\n" + counts.to_csv(sep="\t"),encoding="utf-8")
    pd.DataFrame({"Feature ID":features,"Taxon":[f"d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Genus{i%8};s__Species{i}" for i in range(32)],
                  "Confidence":.99}).to_csv(qiime/"taxonomy.tsv",sep="\t",index=False)
    pd.DataFrame({"SampleID":ids,"pH":rng.normal(6.2,.2,12),"Gas":rng.normal(70,10,12)}).to_csv(qiime/"data.csv",index=False)
    for name,keys in [("KO",[f"K{i:05d}" for i in range(1,21)]),("EC",[f"1.1.1.{i}" for i in range(1,21)])]:
        table = pd.DataFrame(rng.uniform(0,5,size=(32,20)),index=features,columns=keys)
        table.index.name = "ASV"
        table.to_csv(picrust/f"{name}_predicted.tsv.gz",sep="\t",compression="gzip")
    pd.DataFrame({"metadata_NSTI":rng.uniform(.01,.12,32)},index=features).to_csv(picrust/"marker_predicted_and_nsti.tsv.gz",sep="\t",compression="gzip")
    return source


@pytest.mark.parametrize("mode", ["single_factor", "two_factor"])
def test_complete_16s_pipeline(mode, tmp_path):
    source = make_input(tmp_path)
    output = tmp_path / "results"
    stages = SixteenSService.normalize_selected_scripts(PIPELINE_ORDER, mode)
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8", PYTHONUTF8="1")
    for name in stages:
        command = [sys.executable,str(ROOT/"modules/module_16s/scripts"/name),"--input_dir",str(source),
                   "--output_dir",str(output),"--mode",mode]
        if name=="05_picrust2.py": command += ["--app_base_dir",str(ROOT)]
        process = subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True,encoding="utf-8",errors="replace",timeout=240)
        (tmp_path/f"{name}.log").write_text(process.stdout+process.stderr,encoding="utf-8")
        assert process.returncode == 0, (mode,name,process.stdout[-2500:],process.stderr[-2500:])
    assert len(list(output.rglob("*.tsv"))) > 10
    assert len(list(output.rglob("*.png"))) > 5
