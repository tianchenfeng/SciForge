from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
import json
import logging
from pathlib import Path
import shutil
import time
import traceback
from typing import Any, Callable

import numpy as np
import pandas as pd

from modules.module_metabolomics.adapters import GenericMatrixAdapter, WorkflowResultsAdapter
from modules.module_metabolomics.annotation import audit_annotations
from modules.module_metabolomics.correlation import correlate_16s_metabolomics, correlate_features, load_phenotype
from modules.module_metabolomics.enrichment import run_pathway_enrichment
from modules.module_metabolomics.export_results import export_excel, export_html_report, export_manifest, file_sha256, write_text
from modules.module_metabolomics.factorial import run_two_factor
from modules.module_metabolomics.multivariate import Q2_INTERPRETATION_NOTE, Q2_INTERPRETATION_RULES, run_pca, run_plsda
from modules.module_metabolomics.plotting import (
    plot_annotation, plot_correlation_heatmap, plot_filtering, plot_heatmap, plot_intensity,
    plot_interaction_panels, plot_missing, plot_pathway, plot_pca, plot_plsda, plot_qc_rsd,
    plot_scores, plot_top_metabolites, plot_two_factor_heatmap, plot_two_factor_overview,
    plot_correlation_scatterplots, plot_volcano,
)
from modules.module_metabolomics.preprocessing import PreprocessingConfig, preprocess_dataset
from modules.module_metabolomics.qc import compute_qc
from modules.module_metabolomics.statistics import fold_change_rule, significance_class, summarize_significant_features
from modules.module_metabolomics.univariate import run_one_factor
from modules.module_metabolomics.validation import has_blocking_errors, validate_dataset
from modules.module_metabolomics.vendor import load_vendor_oplsda


SCIFORGE_VERSION = "4.0-metabolomics"


@dataclass
class AnalysisConfig:
    input_path: str
    output_parent: str
    analysis_mode: str = "one_factor"
    factor_a: str = "Group"
    factor_b: str = "Factor B"
    metadata_file: str | None = None
    metadata_records: list[dict[str, Any]] | None = None
    annotation_file: str | None = None
    confirmed_input_scale: str = "unknown"
    normalized: bool | None = None
    missing_method: str = "none"
    transformation: str = "none"
    scaling: str = "pareto"
    max_missing_fraction: float = 1.0
    qc_rsd_threshold: float | None = None
    one_factor_method: str = "anova"
    posthoc: str = "tukey"
    fdr_threshold: float = 0.05
    fc_threshold: float = 1.2
    plsda_enabled: bool = True
    plsda_max_components: int = 5
    plsda_variance_threshold: float | None = None
    cv_folds: int = 5
    permutations: int = 200
    random_seed: int = 2026
    heatmap_top_n: int = 30
    phenotype_file: str | None = None
    correlation_method: str = "spearman"
    taxon_table: str | None = None
    title_in_figures: bool = False


def _flatten_figures(value) -> list[str]:
    paths = []
    if isinstance(value, dict):
        for item in value.values(): paths.extend(_flatten_figures(item))
    elif isinstance(value, (list, tuple)):
        for item in value: paths.extend(_flatten_figures(item))
    elif isinstance(value, str): paths.append(value)
    return paths


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(value)).strip("_") or "comparison"


class MetabolomicsAnalysisService:
    @staticmethod
    def detect_adapter(source: str | Path):
        if WorkflowResultsAdapter.can_load(source): return WorkflowResultsAdapter()
        if GenericMatrixAdapter.can_load(source): return GenericMatrixAdapter()
        raise ValueError("Input is neither a supported workflow_results directory nor a generic abundance matrix.")

    @classmethod
    def inspect(cls, source: str | Path, metadata_file: str | None = None, annotation_file: str | None = None):
        adapter = cls.detect_adapter(source)
        dataset = adapter.load(source, metadata_file=metadata_file, annotation_file=annotation_file)
        validation = validate_dataset(dataset)
        return dataset, validation

    @classmethod
    def run(cls, config: AnalysisConfig | dict, progress: Callable[[str, int, str], None] | None = None) -> dict:
        if isinstance(config, dict): config = AnalysisConfig(**config)
        started = time.time(); timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_root = Path(config.output_parent) / f"Metabolomics_Analysis_{timestamp}"
        suffix = 1
        while output_root.exists():
            output_root = Path(config.output_parent) / f"Metabolomics_Analysis_{timestamp}_{suffix}"; suffix += 1
        output_root.mkdir(parents=True)
        log_path = output_root / "analysis.log"
        logging.basicConfig(filename=log_path, level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", force=True, encoding="utf-8")
        warnings: list[str] = []; figures: dict[str, Any] = {}; tables: dict[str, pd.DataFrame] = {}; outputs: dict[str, str] = {}

        def emit(stage, percent, message=""):
            logging.info("[%s %s%%] %s", stage, percent, message)
            if callable(progress): progress(stage, int(percent), message)

        try:
            logging.info("SciForge metabolomics task started")
            emit("Loading data", 2, config.input_path)
            adapter = cls.detect_adapter(config.input_path)
            logging.info("Selected adapter: %s", adapter.name)
            dataset = adapter.load(config.input_path, metadata_file=config.metadata_file, annotation_file=config.annotation_file, scale=config.confirmed_input_scale, normalized=config.normalized)
            if config.metadata_records:
                dataset.samples = pd.DataFrame(config.metadata_records)
            if config.confirmed_input_scale != "unknown": dataset.scale = config.confirmed_input_scale
            if config.analysis_mode in {"one_factor", "two_factor"} and dataset.scale == "unknown":
                raise ValueError("Abundance scale is unknown. Confirm raw, log2, or log10 before formal fold-change analysis.")

            validation = validate_dataset(dataset); tables["Data_Validation"] = validation
            qc_dir = output_root / "01_data_qc"; qc_dir.mkdir()
            validation.to_excel(qc_dir / "data_validation.xlsx", index=False)
            if has_blocking_errors(validation):
                errors = validation.loc[validation["Status"].eq("ERROR"), "Check"].tolist()
                raise ValueError("Blocking data validation errors: " + ", ".join(errors))
            warnings.extend(validation.loc[validation["Status"].eq("WARNING"), "Message"].astype(str).tolist())
            for warning in warnings: logging.warning("%s", warning)

            manifest_dir = output_root / "00_manifest"; manifest_dir.mkdir()
            pd.DataFrame(dataset.freeze_metadata()).to_csv(manifest_dir / "sample_metadata_snapshot.csv", index=False, encoding="utf-8-sig")
            emit("Preprocessing", 8, "Applying explicit preprocessing settings")
            before = dataset.abundance.copy()
            prep_config = PreprocessingConfig(missing_method=config.missing_method, transformation=config.transformation, scaling=config.scaling, max_missing_fraction=config.max_missing_fraction, qc_rsd_threshold=config.qc_rsd_threshold, confirmed_input_scale=config.confirmed_input_scale)
            dataset, filtering = preprocess_dataset(dataset, prep_config); tables["Filtering"] = filtering
            prep_dir = output_root / "02_preprocessing"; prep_dir.mkdir()
            missing_by_feature = before.isna().mean(axis=1).rename("Missing fraction").rename_axis("Feature ID").reset_index()
            figures["QC-01"] = plot_missing(pd.DataFrame({"Sample ID": before.columns, "Missing fraction": before.isna().mean().to_numpy()}), missing_by_feature, prep_dir / "Fig_QC-01_missing")
            figures["QC-02"] = plot_intensity(before, dataset.abundance, prep_dir / "Fig_QC-02_abundance_before_after")
            figures["QC-03"] = plot_filtering(filtering, prep_dir / "Fig_QC-03_filtering")

            emit("QC", 14, "Calculating QC stability metrics")
            qc = compute_qc(dataset); tables["QC"] = qc["qc_table"]; tables["QC_Correlation"] = qc["qc_correlation"]
            tables["QC_Summary"] = pd.DataFrame({"Metric": qc["summary"].keys(), "Value": qc["summary"].values()})
            tables["QC_Summary"].to_csv(qc_dir / "qc_summary.csv", index=False, encoding="utf-8-sig")
            if not qc["qc_table"].empty:
                figures["QC-05"] = plot_correlation_heatmap(qc["qc_correlation"], qc_dir / "Fig_QC-05_correlation", "QC sample correlation")
                figures["QC-06"] = plot_qc_rsd(qc["qc_table"], qc_dir / "Fig_QC-06_RSD")
            else: warnings.append(str(qc["summary"].get("Message")))

            emit("PCA", 20, "Unsupervised structure")
            pca_dir = output_root / "03_pca"; pca_dir.mkdir()
            pca = run_pca(dataset.abundance, dataset.samples, config.scaling)
            pca["scores"].to_csv(pca_dir / "PCA_scores.csv", index=False); pca["loadings"].to_csv(pca_dir / "PCA_loadings.csv", index=False); pca["variance"].to_csv(pca_dir / "PCA_explained_variance.csv", index=False)
            tables.update({"PCA": pca["scores"], "PCA_Loadings": pca["loadings"], "PCA_Variance": pca["variance"]}); figures.update(plot_pca(pca, pca_dir))
            figures["QC-04"] = plot_scores(pca["scores"], pca["variance"], "PC1", "PC2", qc_dir / "Fig_QC-04_PCA_scores")

            pls = None
            if config.plsda_enabled:
                emit("PLS-DA", 25, "Nested cross-validation and permutation testing")
                pls_dir = output_root / "04_plsda"; pls_dir.mkdir()
                pls_input = dataset.modeling_abundance if dataset.modeling_abundance is not None else dataset.abundance
                pls = run_plsda(
                    pls_input, dataset.samples, config.scaling, config.plsda_max_components,
                    config.cv_folds, config.permutations, config.random_seed,
                    progress=lambda value: emit("PLS-DA permutation", 25 + int(value * 15)),
                    missing_method=config.missing_method, transformation=config.transformation,
                    input_scale=dataset.modeling_scale or config.confirmed_input_scale,
                    max_missing_fraction=config.max_missing_fraction,
                    variance_threshold=config.plsda_variance_threshold,
                )
                pls["scores"].to_csv(pls_dir / "PLSDA_scores.csv", index=False); pls["oof_predictions"].to_csv(pls_dir / "PLSDA_oof_predictions.csv", index=False); pls["vip"].to_csv(pls_dir / "PLSDA_VIP.csv", index=False); pls["cv"].to_csv(pls_dir / "PLSDA_cv.csv", index=False); pls["permutation"].to_csv(pls_dir / "PLSDA_permutation.csv", index=False)
                tables.update({"PLSDA": pls["metrics"], "PLSDA_OOF": pls["oof_predictions"], "PLSDA_VIP": pls["vip"], "PLSDA_CV": pls["cv"], "PLSDA_Permutation": pls["permutation"]}); figures.update(plot_plsda(pls, pls_dir))

            annotated, annotation_stats = audit_annotations(dataset.features); tables["Feature_Metadata"] = annotated; tables["Annotations"] = annotated[annotated["Identified"]].copy()
            anno_dir = output_root / "09_annotation"; anno_dir.mkdir(); figures.update(plot_annotation(annotation_stats, anno_dir)); annotated.to_csv(anno_dir / "annotation_audit.csv", index=False, encoding="utf-8-sig")
            feature_info = annotated.set_index("Feature ID"); labels = feature_info["Metabolite name"].replace("", pd.NA).fillna(feature_info.index.to_series()).to_dict()
            heat_dir = output_root / "08_heatmap"; heat_dir.mkdir()
            annotated_ids = annotated.loc[annotated["Identified"], "Feature ID"].astype(str)
            biological_ids = dataset.biological_samples
            ranked_annotated = dataset.abundance.loc[dataset.abundance.index.intersection(annotated_ids), biological_ids].var(axis=1).sort_values(ascending=False).index.astype(str).tolist()
            if ranked_annotated:
                figures["HEAT-02"] = plot_heatmap(dataset.abundance, dataset.samples, ranked_annotated, labels, heat_dir / "Fig_HEAT-02_annotated_metabolites", config.heatmap_top_n)

            overall = contrasts = two_factor = simple = pd.DataFrame(); differential_frames = []; pathway_frames = []
            if config.analysis_mode == "one_factor":
                emit("Univariate", 42, "Overall tests and post-hoc contrasts")
                uni_dir = output_root / "05_univariate"; uni_dir.mkdir()
                overall, contrasts = run_one_factor(dataset.abundance, dataset.samples, config.factor_a, config.one_factor_method, config.posthoc, dataset.scale)
                tables["Overall_Test"] = overall; tables["Contrasts"] = contrasts
                overall.to_csv(uni_dir / "overall_test.csv", index=False); contrasts.to_csv(uni_dir / "contrasts.csv", index=False)
                diff_dir = output_root / "07_differential"; diff_dir.mkdir(); pathway_dir = output_root / "10_pathway"; pathway_dir.mkdir()
                for contrast_name, block in contrasts.groupby("Contrast"):
                    block = block.copy(); block["Significance class"] = significance_class(block["FDR"], block["FC"], config.fdr_threshold, config.fc_threshold)
                    merged = block.merge(annotated, on="Feature ID", how="left")
                    fc_note = fold_change_rule(dataset.scale, config.fc_threshold, dataset.provenance.get("preprocessing", {}).get("transformation_pseudocount"))
                    merged["FC direction"] = fc_note["contrast_direction"]
                    merged["FC formula"] = fc_note["fc_formula"]
                    merged["log2FC formula"] = fc_note["log2fc_formula"]
                    differential_frames.append(merged)
                    safe = _safe_name(contrast_name); merged.to_csv(diff_dir / f"differential_{safe}.csv", index=False, encoding="utf-8-sig")
                    figures[f"DIFF-FDR-{safe}"] = plot_volcano(merged, diff_dir / f"Fig_DIFF-01_{safe}_FDR", "FDR", False, config.fdr_threshold, config.fc_threshold)
                    figures[f"DIFF-P-{safe}"] = plot_volcano(merged, diff_dir / f"Fig_DIFF-02_{safe}_exploratory_rawP", "P", True, .05, config.fc_threshold)
                    figures[f"DIFF-TOP-{safe}"] = plot_top_metabolites(dataset.abundance, dataset.samples, merged, labels, diff_dir / f"Fig_DIFF-04_{safe}_top_metabolites")
                    significant_ids = merged.loc[merged["Significance class"].isin(["Up", "Down"]), "Feature ID"].astype(str).tolist()
                    if significant_ids: figures[f"HEAT-{safe}"] = plot_heatmap(dataset.abundance, dataset.samples, significant_ids, labels, heat_dir / f"Fig_HEAT-01_{safe}", config.heatmap_top_n)
                    enrichment = run_pathway_enrichment(significant_ids, annotated)
                    if not enrichment.empty: enrichment.insert(0, "Contrast", contrast_name); pathway_frames.append(enrichment); enrichment.to_csv(pathway_dir / f"pathway_{safe}.csv", index=False, encoding="utf-8-sig")
                    if not enrichment.empty and (enrichment["FDR"] < config.fdr_threshold).any(): figures[f"PATH-{safe}"] = plot_pathway(enrichment, pathway_dir / f"Fig_PATH-01_{safe}")
                    else: warnings.append(f"{contrast_name}: No pathways remained significant after FDR correction.")
                if differential_frames:
                    membership = pd.concat(differential_frames, ignore_index=True)
                    membership = membership[membership["Significance class"].isin(["Up", "Down"])]
                    if not membership.empty:
                        comparison_names = list(map(str, pd.unique(membership["Contrast"])))
                        rows = []
                        for feature_id, block in membership.groupby("Feature ID"):
                            present = set(block["Contrast"].astype(str)); rows.append({"Feature ID": feature_id, "Intersection": " & ".join(name for name in comparison_names if name in present), **{name: name in present for name in comparison_names}})
                        intersections = pd.DataFrame(rows); intersections.to_csv(diff_dir / "differential_intersection_table.csv", index=False, encoding="utf-8-sig"); tables["Intersections"] = intersections
            else:
                emit("Two-factor", 42, "A, B, interaction and model contrasts")
                two_dir = output_root / "06_two_factor"; two_dir.mkdir()
                two_factor, simple = run_two_factor(dataset.abundance, dataset.samples, config.factor_a, config.factor_b, config.fdr_threshold, progress=lambda value: emit("Two-factor", 42 + int(value * 20)))
                tables["TwoFactor_A"] = two_factor[[c for c in two_factor if c in ["Feature ID", "F_A", "P_A", "q_A", "effect_A", "Model"]]]
                tables["TwoFactor_B"] = two_factor[[c for c in two_factor if c in ["Feature ID", "F_B", "P_B", "q_B", "effect_B", "Model"]]]
                tables["TwoFactor_Interaction"] = two_factor[[c for c in two_factor if c in ["Feature ID", "F_AxB", "P_AxB", "q_AxB", "effect_AxB", "Model"]]]; tables["SimpleEffects"] = simple
                two_factor.to_csv(two_dir / "two_factor_results.csv", index=False); simple.to_csv(two_dir / "simple_effects.csv", index=False)
                for effect, number in [("A", "01"), ("B", "02"), ("AxB", "03")]: figures[f"2F-{effect}"] = plot_two_factor_overview(two_factor, effect, two_dir / f"Fig_2F-{number}_{effect}")
                figures["2F-interaction"] = plot_interaction_panels(dataset.abundance, dataset.samples, two_factor, config.factor_a, config.factor_b, labels, two_dir / "Fig_2F-04_interaction_top", fdr_threshold=config.fdr_threshold)
                figures["2F-heatmap"] = plot_two_factor_heatmap(two_factor, labels, two_dir / "Fig_2F-05_effect_heatmap")

            differential = pd.concat(differential_frames, ignore_index=True) if differential_frames else pd.DataFrame(); pathways = pd.concat(pathway_frames, ignore_index=True) if pathway_frames else pd.DataFrame()
            tables["Differential"] = differential; tables["Pathways"] = pathways

            correlations = pd.DataFrame(); multiomics = pd.DataFrame()
            selected_features = differential.loc[differential.get("Significance class", pd.Series(dtype=object)).isin(["Up", "Down"]), "Feature ID"].drop_duplicates().astype(str).tolist() if not differential.empty else None
            if config.phenotype_file:
                emit("Correlation", 66, "Metabolite × phenotype")
                corr_dir = output_root / "11_correlation"; corr_dir.mkdir(); phenotype = load_phenotype(config.phenotype_file)
                correlations, corr_matrix, q_matrix = correlate_features(dataset.abundance, phenotype, config.correlation_method, selected_features)
                correlations.to_csv(corr_dir / "correlation_long.csv", index=False); correlations.assign(Significant=correlations["FDR"] < config.fdr_threshold).to_csv(corr_dir / "network_ready.csv", index=False)
                tables["Correlations"] = correlations; figures["COR-01"] = plot_correlation_heatmap(corr_matrix, corr_dir / "Fig_COR-01_heatmap", "Metabolite × phenotype")
                figures["COR-02"] = plot_correlation_scatterplots(dataset.abundance, phenotype, correlations, labels, corr_dir / "Fig_COR-02_top_scatterplots")
            if config.taxon_table:
                from modules.module_metabolomics.adapters.base import read_table_smart
                multi_dir = output_root / "12_multiomics"; multi_dir.mkdir(); taxon = read_table_smart(config.taxon_table, low_memory=False); taxon = taxon.set_index(taxon.columns[0])
                multiomics = correlate_16s_metabolomics(taxon, dataset.abundance)
                multiomics.to_csv(multi_dir / "correlation_long.csv", index=False); multiomics.to_csv(multi_dir / "network_ready.csv", index=False)
                multi_matrix = multiomics.pivot(index="Taxon", columns="Feature ID", values="Spearman")
                top_rows = multiomics.dropna(subset=["FDR"]).sort_values("FDR").head(1000)
                if not top_rows.empty:
                    multi_matrix = multi_matrix.loc[multi_matrix.index.intersection(top_rows["Taxon"]), multi_matrix.columns.intersection(top_rows["Feature ID"])]
                    figures["MULTI-COR"] = plot_correlation_heatmap(multi_matrix, multi_dir / "Fig_MULTI-01_correlation_heatmap", "16S × metabolomics Spearman")

            vendor_opls = load_vendor_oplsda(dataset)
            vendor_sheet_names = {"model": "Vendor_OPLS_Model", "scores": "Vendor_OPLS_Scores", "vip": "Vendor_OPLS_VIP", "permutation": "Vendor_OPLS_Permutation"}
            for key, table in vendor_opls.items():
                tables[vendor_sheet_names.get(key, f"Vendor_OPLS_{key}")] = table

            emit("Export", 78, "Writing reproducible result package")
            significant = summarize_significant_features(differential)
            fc_config = fold_change_rule(dataset.scale, config.fc_threshold, dataset.provenance.get("preprocessing", {}).get("transformation_pseudocount"))
            pls_metrics = None if pls is None else pls["metrics"].iloc[0].to_dict()
            q2_statement = "PLS-DA was not run."
            if pls_metrics is not None:
                q2_statement = (
                    f"PLS-DA 的 out-of-fold Q² 为 {float(pls_metrics['Q2']):.4f}，表明{pls_metrics['Predictive interpretation']}。"
                    f"置换检验状态为 {pls_metrics['Permutation status']}；置换检验通过仅表示模型优于随机标签，不等于具有良好的预测性能。"
                )
            summary = {
                "Biological samples": len(dataset.biological_samples), "QC samples": len(dataset.qc_samples), "Features": len(dataset.abundance),
                "Annotated metabolites": int(annotated["Identified"].sum()), "KEGG mapped": int((~annotated["KEGG"].fillna("").astype(str).isin(["", "-", "nan"])).sum()),
                "QC RSD <30%": qc["summary"].get("RSD <30%", "Unavailable"), "PCA PC1": float(pca["variance"].iloc[0]["Explained variance"]), "PCA PC2": float(pca["variance"].iloc[1]["Explained variance"]),
                "Significant feature-contrast records": significant["records"], "Unique significant features": significant["unique_features"], "Up feature-contrast records": significant["up"], "Down feature-contrast records": significant["down"], "Significant pathways": int((pathways.get("FDR", pd.Series(dtype=float)) < config.fdr_threshold).sum()),
                "Two-factor A hits": int((two_factor.get("q_A", pd.Series(dtype=float)) < config.fdr_threshold).sum()), "Two-factor B hits": int((two_factor.get("q_B", pd.Series(dtype=float)) < config.fdr_threshold).sum()), "Two-factor interaction hits": int((two_factor.get("q_AxB", pd.Series(dtype=float)) < config.fdr_threshold).sum()),
                "PLS-DA OOF Q2": None if pls_metrics is None else float(pls_metrics["Q2"]), "PLS-DA predictive strength": "not run" if pls_metrics is None else pls_metrics["Predictive strength"], "PLS-DA interpretation": q2_statement, "PLS-DA permutation status": "Not run" if pls_metrics is None else pls_metrics["Permutation status"],
            }
            tables["Summary"] = pd.DataFrame({"Metric": summary.keys(), "Value": summary.values()}); tables["Samples"] = dataset.samples
            q2_rules = {"scope": Q2_INTERPRETATION_NOTE, "rules": Q2_INTERPRETATION_RULES}
            parameter_values = {**asdict(config), "fold_change": fc_config, "Q2 interpretation rules": q2_rules}
            tables["Parameters"] = pd.DataFrame({"Parameter": parameter_values.keys(), "Value": [json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for v in parameter_values.values()]})
            table_dir = output_root / "tables"; table_dir.mkdir(); excel = export_excel(tables, table_dir / "metabolomics_analysis.xlsx"); outputs["excel"] = excel

            fact_lines = [f"- {key}: {value}" for key, value in summary.items()]
            if not differential.empty:
                for name, block in differential.groupby("Contrast"):
                    count = int(block["Significance class"].isin(["Up", "Down"]).sum()); annotated_count = int(block.loc[block["Significance class"].isin(["Up", "Down"]), "Identified"].fillna(False).sum())
                    fact_lines.append(f"- {name}: {count} feature-contrast records met FDR<{config.fdr_threshold} and FC threshold {config.fc_threshold}; {annotated_count} were annotated.")
            report_dir = output_root / "report"; report_dir.mkdir(); outputs["summary_md"] = write_text(report_dir / "analysis_summary.md", "# Analysis Summary\n\n" + "\n".join(fact_lines) + "\n")
            normalization_text = "existing normalized abundance" if dataset.normalized is True else "no additional normalization"
            methods = [f"Input data were imported with the {adapter.name} adapter.", f"Normalization: {normalization_text}; missing values: {config.missing_method}; transformation: {config.transformation}; scaling: {config.scaling}.", f"Multiple testing was controlled using Benjamini-Hochberg FDR; the formal threshold was q<{config.fdr_threshold}."]
            if config.analysis_mode == "one_factor": methods.append(f"One-factor analysis used {config.one_factor_method} with {config.posthoc} post-hoc comparisons and an FC threshold of {config.fc_threshold}. Contrasts are labeled TARGET vs REFERENCE; {fc_config['fc_formula']}; {fc_config['log2fc_formula']}; FC > 1 means TARGET is higher than REFERENCE. Zero handling: {fc_config['zero_handling']}")
            else: methods.append(f"Two-factor analysis modeled {config.factor_a} * {config.factor_b}; balanced designs used Type-II and unbalanced designs used sum-contrast Type-III ANOVA. Significant interactions were followed by model-based simple effects.")
            methods.append("PCA was unsupervised. " + (f"PLS-DA used strict nested stratified CV with fold-local missing-fraction/variance feature selection, imputation, transformation, scaling, and {config.permutations} full nested label permutations. Q² used pooled outer-test predictions. {q2_statement} The predictive-strength bands are a SciForge conservative reporting rule, not a universal metabolomics industry standard. VIP was interpretive and not a significance test." if pls is not None else "PLS-DA was not run."))
            if not pathways.empty: methods.append("Pathway over-representation used explicit KEGG compound identifiers and hypergeometric tests with BH-FDR correction.")
            if config.phenotype_file: methods.append(f"Metabolite–phenotype associations used {config.correlation_method} correlation with BH-FDR correction across tested pairs.")
            if config.taxon_table: methods.append("The lightweight 16S–metabolomics analysis used Spearman correlation with BH-FDR correction across taxon–feature pairs.")
            outputs["methods"] = write_text(report_dir / "methods_metabolomics.txt", "\n".join(methods) + "\n")
            pngs = sorted({p for p in _flatten_figures(figures) if p.endswith(".png")})
            legends = [f"- {Path(path).stem}. Generated from the current analysis using n={len(dataset.biological_samples)} biological samples; scaling={config.scaling}; formal significance=q<{config.fdr_threshold} and FC threshold={config.fc_threshold} where applicable." for path in pngs]
            outputs["legends"] = write_text(report_dir / "figure_legends.md", "# Figure Legends\n\n" + "\n".join(legends) + "\n")
            sections = [("Project summary", f"Adapter={adapter.name}; design={config.analysis_mode}."), ("Data quality", f"{len(dataset.abundance)} features and {len(dataset.qc_samples)} QC samples were analyzed."), ("Preprocessing", methods[1]), ("PCA", f"PC1 and PC2 explained {summary['PCA PC1']:.1%} and {summary['PCA PC2']:.1%}."), ("PLS-DA", q2_statement), ("Univariate", f"Run with {config.one_factor_method}." if config.analysis_mode == "one_factor" else "Not run."), ("Fold change", f"Direction: TARGET vs REFERENCE. FC formula: {fc_config['fc_formula']}. log2FC formula: {fc_config['log2fc_formula']}."), ("Two-factor", f"Modeled {config.factor_a} × {config.factor_b}." if config.analysis_mode == "two_factor" else "Not run."), ("Differential metabolites", f"{significant['records']} significant feature-contrast records involving {significant['unique_features']} unique features."), ("Annotation", f"{summary['Annotated metabolites']} features had annotation confidence metadata."), ("Pathways", f"FDR-significant pathways: {summary['Significant pathways']}."), ("Correlation", "Run" if config.phenotype_file else "Not run"), ("Methods", " ".join(methods))]
            outputs["html"] = export_html_report(summary, sections, pngs, warnings, report_dir / "metabolomics_report.html")

            paper_dir = output_root / "figures"; paper_dir.mkdir()
            for index, path in enumerate(pngs, 1): shutil.copy2(path, paper_dir / f"Fig{index:02d}_{Path(path).name}")
            input_hashes = {}
            for key, value in dataset.provenance.items():
                if isinstance(value, str) and Path(value).is_file(): input_hashes[key] = file_sha256(value)
            pls_manifest = None if pls_metrics is None else {**pls_metrics, "summary": q2_statement, "Q2 interpretation rules": q2_rules, "VIP role": "interpretive only; not a formal significance criterion"}
            manifest = {"SciForge version": SCIFORGE_VERSION, "timestamp": timestamp, "runtime_seconds": time.time() - started, "input paths": dataset.provenance, "input hashes": input_hashes, "adapter": adapter.name, "sample metadata": dataset.freeze_metadata(), "groups": dataset.samples["Group"].value_counts().to_dict(), "factor A": config.factor_a, "factor B": config.factor_b if config.analysis_mode == "two_factor" else None, "preprocessing": dataset.provenance.get("preprocessing", prep_config.to_dict()), "statistics": {"mode": config.analysis_mode, "one_factor_method": config.one_factor_method, "posthoc": config.posthoc, "two_factor": "Type-II if balanced; sum contrasts + Type-III if unbalanced"}, "FDR method": "Benjamini-Hochberg", "thresholds": {"FDR": config.fdr_threshold, "FC": config.fc_threshold}, "fold_change": fc_config, "analysis summary": summary, "significant features": significant, "PCA scaling": config.scaling, "PLSDA": pls_manifest, "random seed": config.random_seed, "output paths": outputs, "warnings": sorted(set(warnings))}
            outputs["manifest"] = export_manifest(manifest, manifest_dir / "analysis_manifest.json")
            logging.info("Analysis summary: %s", json.dumps(summary, ensure_ascii=False))
            emit("Complete", 100, str(output_root)); logging.info("Final status: success; runtime %.2fs", time.time() - started)
            return {"success": True, "output_dir": str(output_root), "summary": summary, "outputs": outputs, "warnings": sorted(set(warnings)), "runtime_seconds": time.time() - started}
        except Exception as exc:
            logging.error("Analysis failed: %s\n%s", exc, traceback.format_exc()); emit("Failed", 100, str(exc))
            return {"success": False, "output_dir": str(output_root), "message": str(exc), "log": str(log_path), "runtime_seconds": time.time() - started}
