from copy import deepcopy
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from modules.module_feed_formula.service import FeedFormulaService
from modules.module_feed_formula.validators import FeedFormulaValidator
from modules.module_metabolomics.data_model import MetabolomicsDataset
from modules.module_metabolomics.service import AnalysisConfig, MetabolomicsAnalysisService
from modules.module_metabolomics.validation import validate_dataset
from ui.page_feed_formula import PageFeedFormula


@pytest.mark.parametrize("module", ["gas_fit", "silage_eval", "rnaseq_analysis"])
def test_original_script_acceptance(module, monkeypatch, tmp_path):
    script = importlib.import_module(f"modules.module_{module}.tests.test_{module}")
    monkeypatch.setattr(script, "__file__", str(tmp_path / "acceptance.py"))
    script.main()


def ingredients_and_constraints():
    return deepcopy(PageFeedFormula.DEFAULT_ROWS), {"cp": {"enable":True,"target":17.5,"error":1.},
        "starch":{"enable":True,"target":24.,"error":2.}, "ratio":{"enable":True,"target":50.,"error":5.}}


def test_feed_valid_solution_and_impossible_constraints():
    ingredients, constraints = ingredients_and_constraints()
    result = FeedFormulaService.run_formula(ingredients, constraints)
    assert result["success"], result["message"]
    assert result["mix"] and result["stats"]
    constraints["cp"].update(target=500.,error=0.)
    impossible = FeedFormulaService.run_formula(ingredients, constraints)
    assert impossible["success"] is False and impossible["message"]


@pytest.mark.parametrize("field", FeedFormulaValidator.NUMERIC_FIELDS)
@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf"), -1., "not numeric"])
def test_feed_invalid_ingredient_numbers(field, value):
    ingredients, constraints = ingredients_and_constraints()
    ingredients[0][field] = value
    with pytest.raises(ValueError): FeedFormulaValidator.validate(ingredients, constraints)


@pytest.mark.parametrize("field", ["target", "error"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1., "not numeric"])
def test_feed_invalid_constraint_numbers(field, value):
    ingredients, constraints = ingredients_and_constraints()
    constraints["cp"][field] = value
    with pytest.raises(ValueError): FeedFormulaValidator.validate(ingredients, constraints)


@pytest.mark.parametrize("case", ["empty", "duplicate", "missing_name", "invalid_type", "no_enabled", "no_roughage"])
def test_feed_structure_validation(case):
    ingredients, constraints = ingredients_and_constraints()
    if case == "empty": ingredients = []
    elif case == "duplicate": ingredients[1]["name"] = ingredients[0]["name"]
    elif case == "missing_name": ingredients[0]["name"] = ""
    elif case == "invalid_type": ingredients[0]["type"] = "unknown"
    elif case == "no_enabled":
        for value in constraints.values(): value["enable"] = False
    elif case == "no_roughage":
        for value in ingredients: value["type"] = "concentrate"
    with pytest.raises(ValueError):
        FeedFormulaValidator.validate(ingredients, constraints)
        FeedFormulaValidator.validate_roughage_presence(ingredients, constraints)


def make_metabolomics(tmp_path):
    rng = np.random.default_rng(2626)
    ids = [f"S{i}" for i in range(16)]
    matrix = rng.lognormal(4, .3, (20,16))
    matrix[:4,8:] *= 2
    abundance = pd.DataFrame(matrix,index=[f"F{i}" for i in range(20)],columns=ids)
    metadata = pd.DataFrame([{"Sample ID":s,"Sample type":"biological", "Group":"C" if i<8 else "T",
         "Factor A":"A1" if i<8 else "A2", "Factor B":"B1" if (i//4)%2==0 else "B2"} for i,s in enumerate(ids)])
    features = pd.DataFrame({"Feature ID":abundance.index})
    return MetabolomicsDataset(abundance,metadata,features,scale="raw")


@pytest.mark.parametrize("case", ["infinity","huge","tiny","duplicate_feature","duplicate_sample","missing_group"])
def test_metabolomics_invalid_data_is_blocked(case, tmp_path):
    dataset = make_metabolomics(tmp_path)
    if case == "infinity": dataset.abundance.iloc[0,0] = np.inf
    elif case == "huge": dataset.abundance.iloc[0,0] = 1e300
    elif case == "tiny": dataset.abundance.iloc[0,0] = 1e-300
    elif case == "duplicate_feature": dataset.abundance.index = ["same"]*len(dataset.abundance)
    elif case == "duplicate_sample": dataset.samples.loc[0,"Sample ID"] = dataset.samples.loc[1,"Sample ID"]
    elif case == "missing_group": dataset.samples.loc[0,"Group"] = ""
    checks = validate_dataset(dataset)
    assert checks["Status"].eq("ERROR").any()


@pytest.mark.parametrize("mode", ["one_factor", "two_factor"])
def test_metabolomics_end_to_end_both_designs(mode, tmp_path):
    dataset = make_metabolomics(tmp_path)
    frame = dataset.abundance.rename_axis("Feature ID").reset_index()
    source = tmp_path / "matrix.csv"
    frame.to_csv(source,index=False)
    result = MetabolomicsAnalysisService.run(AnalysisConfig(
        input_path=str(source),output_parent=str(tmp_path),metadata_records=dataset.samples.to_dict("records"),
        analysis_mode=mode,factor_a="Factor A",factor_b="Factor B",confirmed_input_scale="raw",
        permutations=2, cv_folds=3,plsda_max_components=2, heatmap_top_n=10))
    assert result["success"], result.get("message")
    output = Path(result["output_dir"])
    assert (output/"tables/metabolomics_analysis.xlsx").exists()
    assert (output/"report/metabolomics_report.html").exists()
    manifest = json.loads((output/"00_manifest/analysis_manifest.json").read_text(encoding="utf-8"))
    assert manifest["statistics"]["mode"] == mode


def test_metabolomics_invalid_scale_fails_cleanly(tmp_path):
    dataset = make_metabolomics(tmp_path)
    source = tmp_path/"matrix.csv"
    dataset.abundance.rename_axis("Feature ID").reset_index().to_csv(source,index=False)
    result = MetabolomicsAnalysisService.run(AnalysisConfig(input_path=str(source),output_parent=str(tmp_path)))
    assert result["success"] is False
    assert "scale" in result["message"].lower()
