"""Generated stress inputs exercise the public service and the real document exporters."""
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest
from scipy.stats import f_oneway

from modules.module_stats.service import run_one_way_batch_analysis, run_two_way_batch_analysis
from modules.module_stats.stats_core import run_one_way_anova, run_two_way_anova


def frame_for(mode, seed=51):
    rng = np.random.default_rng(seed)
    rows = [{"Group": a, "Time": b, "value": 10 + ai*2 + bi + rng.normal(0,.3)}
            for ai,a in enumerate(["Control", "Treatment"]) for bi,b in enumerate(["Day 1", "Day 2"])
            for _ in range(4)]
    frame = pd.DataFrame(rows)
    return frame.drop(columns="Time") if mode == "one" else frame


def run(frame, mode, tmp_path, extension="csv"):
    source = tmp_path / f"合成测试数据.{extension}"
    if extension == "xlsx": frame.to_excel(source, index=False)
    else: frame.to_csv(source, index=False, sep="\t" if extension in {"tsv", "txt"} else ",")
    service = run_one_way_batch_analysis if mode == "one" else run_two_way_batch_analysis
    return service(source, tmp_path / "结果目录")


def assert_package(result):
    assert result["success"], result["message"]
    assert result["metrics"]
    for key in ("word_path", "excel_path"):
        with ZipFile(result[key]) as archive:
            assert archive.testzip() is None
    table = pd.read_excel(result["excel_path"], sheet_name="summary_long")
    assert np.isfinite(pd.to_numeric(table["mean"]).to_numpy(dtype=float)).all()
    text = result["paper_table"].astype(str).to_string().lower()
    assert "nan" not in text and "inf" not in text


@pytest.mark.parametrize("mode", ["one", "two"])
@pytest.mark.parametrize("extension", ["csv", "tsv", "txt", "xlsx"])
def test_all_input_formats(mode, extension, tmp_path):
    assert_package(run(frame_for(mode), mode, tmp_path, extension))


@pytest.mark.parametrize("mode", ["one", "two"])
@pytest.mark.parametrize("case", ["zero", "constant", "negative", "large", "small", "mixed_text", "missing_values",
                                    "unbalanced", "unicode", "symbols", "long_names", "numeric_strings", "many_metrics"])
def test_supported_boundaries(mode, case, tmp_path):
    frame = frame_for(mode)
    if case == "zero": frame["value"] = 0.
    elif case == "constant": frame["value"] = 7.
    elif case == "negative": frame["value"] *= -1
    elif case == "large": frame["value"] *= 1e100
    elif case == "small": frame["value"] *= 1e-100
    elif case == "mixed_text":
        frame["value"] = frame["value"].astype(object)
        frame.loc[0, "value"] = "below detection"
        frame["备注"] = "not a metric"
    elif case == "missing_values": frame.loc[[0, 5], "value"] = np.nan
    elif case == "unbalanced": frame = frame.drop(index=[0,1,5])
    elif case == "unicode": frame["Group"] = frame["Group"].map({"Control":"对照组", "Treatment":"处理组"})
    elif case == "symbols": frame = frame.rename(columns={"Group": '因素 A "种类"', "Time": "Time / 天", "value": "NH3-N (mg/L)%"})
    elif case == "long_names": frame = frame.rename(columns={"value": "实验指标_" + "very_long_name_"*12})
    elif case == "numeric_strings": frame["value"] = frame["value"].map(str)
    elif case == "many_metrics":
        for i in range(12): frame[f"metric_{i}"] = frame["value"] + i*.2
    assert_package(run(frame, mode, tmp_path))


@pytest.mark.parametrize("mode", ["one", "two"])
@pytest.mark.parametrize("case", ["infinity", "negative_infinity", "huge", "tiny", "all_nan", "all_text", "missing_group",
                                    "blank_group", "single_group", "empty", "no_metrics", "duplicate_clean_names"])
def test_invalid_inputs_return_actionable_failure(mode, case, tmp_path):
    frame = frame_for(mode)
    if case == "infinity": frame.loc[0,"value"] = np.inf
    elif case == "negative_infinity": frame.loc[0,"value"] = -np.inf
    elif case == "huge": frame["value"] *= 1e290
    elif case == "tiny": frame["value"] *= 1e-290
    elif case == "all_nan": frame["value"] = np.nan
    elif case == "all_text": frame["value"] = "bad"
    elif case == "missing_group": frame.loc[0,"Group"] = None
    elif case == "blank_group": frame.loc[0,"Group"] = "   "
    elif case == "single_group": frame["Group"] = "only"
    elif case == "empty": frame = frame.iloc[:0]
    elif case == "no_metrics": frame = frame.drop(columns="value")
    elif case == "duplicate_clean_names": frame[" value "] = frame["value"]
    result = run(frame, mode, tmp_path)
    assert result["success"] is False
    assert result["message"].strip()
    assert not result["word_path"] and not result["excel_path"]


@pytest.mark.parametrize("mode", ["one", "two"])
@pytest.mark.parametrize("case", ["missing", "empty", "corrupt_excel", "unsupported"])
def test_file_failures(mode, case, tmp_path):
    source = tmp_path / ("input.xlsx" if case == "corrupt_excel" else "input.exe" if case == "unsupported" else "input.csv")
    if case != "missing": source.write_bytes(b"not an excel file" if case == "corrupt_excel" else b"")
    service = run_one_way_batch_analysis if mode == "one" else run_two_way_batch_analysis
    result = service(source, tmp_path / "out")
    assert result["success"] is False and result["message"]


@pytest.mark.parametrize("seed", range(10))
def test_one_way_against_independent_scipy_reference(seed):
    frame = frame_for("one", seed)
    result = run_one_way_anova(frame, "Group", "value")
    reference = f_oneway(*(g["value"].to_numpy() for _,g in frame.groupby("Group")))
    assert result["anova_table"].loc["C(Group)", "F"] == pytest.approx(reference.statistic, rel=1e-10)
    assert result["anova_table"].loc["C(Group)", "PR(>F)"] == pytest.approx(reference.pvalue, rel=1e-10)


@pytest.mark.parametrize("mode", ["one", "two"])
@pytest.mark.parametrize("scale,offset", [(1.,100.), (1000.,0.), (.001,0.), (-3.,2.)])
def test_anova_p_values_invariant_under_unit_conversion(mode, scale, offset):
    frame = frame_for(mode)
    analyze = (lambda f: run_one_way_anova(f,"Group","value")) if mode == "one" else (lambda f: run_two_way_anova(f,"Group","Time","value"))
    first = analyze(frame)["anova_table"]["PR(>F)"].dropna()
    frame["value"] = frame["value"]*scale + offset
    second = analyze(frame)["anova_table"]["PR(>F)"].dropna()
    np.testing.assert_allclose(first, second, rtol=1e-8, atol=1e-14)
