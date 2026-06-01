from __future__ import annotations

from pathlib import Path

import pandas as pd

from modules.module_rnaseq_analysis.annotation_merge import merge_annotation
from modules.module_rnaseq_analysis.candidate_ranker import rank_candidates
from modules.module_rnaseq_analysis.deg_filter import filter_deg, split_up_down
from modules.module_rnaseq_analysis.enrichment import standardize_enrichment_table
from modules.module_rnaseq_analysis.export_results import export_excel
from modules.module_rnaseq_analysis.gene_set_ops import parse_keywords
from modules.module_rnaseq_analysis.pathway_gene_map import build_pathway_gene_map
from modules.module_rnaseq_analysis.plotting import draw_candidate_heatmap, draw_enrichment_barplot, draw_volcano_plot
from modules.module_rnaseq_analysis.preprocess import (
    load_annotation,
    load_expression_matrix,
    read_deg_table,
    standardize_deg_table,
)
from modules.module_rnaseq_analysis.word_export import export_word_report


class RNASeqAnalysisService:
    @staticmethod
    def run(
        deg_file: str,
        output_dir: str,
        log2fc_threshold: float = 1.0,
        padj_threshold: float = 0.05,
        keywords: str | list[str] | None = None,
        annotation_file: str | None = None,
        go_file: str | None = None,
        kegg_file: str | None = None,
        expression_file: str | None = None,
        logger=None,
    ) -> dict:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        def log(message: str):
            if callable(logger):
                logger(message)

        keyword_list = parse_keywords(keywords) if isinstance(keywords, str) else (keywords or [])
        log(f"Reading DEG table: {deg_file}")
        raw = read_deg_table(deg_file)
        all_deg, detected = standardize_deg_table(raw)
        filtered_deg = filter_deg(all_deg, log2fc_threshold, padj_threshold)
        up_genes, down_genes = split_up_down(filtered_deg)

        annotation = load_annotation(annotation_file)
        deg_with_annotation = merge_annotation(all_deg, annotation)
        filtered_with_annotation = merge_annotation(filtered_deg, annotation)
        up_with_annotation = merge_annotation(up_genes, annotation)
        down_with_annotation = merge_annotation(down_genes, annotation)

        go_sets = standardize_enrichment_table(go_file, "GO", keyword_list, padj_threshold)
        kegg_sets = standardize_enrichment_table(kegg_file, "KEGG", keyword_list, padj_threshold)
        pathway_gene_map = build_pathway_gene_map(go_sets["all"], kegg_sets["all"])
        candidates = rank_candidates(filtered_with_annotation, pathway_gene_map, keyword_list)

        volcano = draw_volcano_plot(all_deg, output_path / "volcano_plot.png", log2fc_threshold, padj_threshold)
        top_go_barplot = draw_enrichment_barplot(go_sets["top"], output_path / "top_GO_barplot.png", "Top GO enrichment")
        top_kegg_barplot = draw_enrichment_barplot(kegg_sets["top"], output_path / "top_KEGG_barplot.png", "Top KEGG enrichment")
        expression = load_expression_matrix(expression_file)
        heatmap = draw_candidate_heatmap(expression, candidates, output_path / "candidate_gene_heatmap.png")

        tables = {
            "all_deg": all_deg,
            "filtered_deg": filtered_deg,
            "up_genes": up_genes,
            "down_genes": down_genes,
            "deg_with_annotation": deg_with_annotation,
            "filtered_deg_with_annotation": filtered_with_annotation,
            "up_genes_with_annotation": up_with_annotation,
            "down_genes_with_annotation": down_with_annotation,
            "candidate_genes_ranked": candidates,
            "top_GO": go_sets["top"],
            "significant_GO": go_sets["significant"],
            "keyword_GO": go_sets["keyword"],
            "top_KEGG": kegg_sets["top"],
            "significant_KEGG": kegg_sets["significant"],
            "keyword_KEGG": kegg_sets["keyword"],
            "pathway_gene_map": pathway_gene_map,
        }
        excel = export_excel(tables, output_path / "rnaseq_analysis_results.xlsx")
        summary = {
            "filtered_count": len(filtered_deg),
            "up_count": len(up_genes),
            "down_count": len(down_genes),
        }
        word = export_word_report(summary, candidates, go_sets["top"], kegg_sets["top"], output_path / "rnaseq_analysis_report.docx")
        log(f"Filtered DEG: {len(filtered_deg)}; Up: {len(up_genes)}; Down: {len(down_genes)}")

        return {
            "success": True,
            "message": "RNASeq analysis finished",
            "excel": excel,
            "word": word,
            "volcano": volcano,
            "candidate_heatmap": heatmap,
            "top_go_barplot": top_go_barplot,
            "top_kegg_barplot": top_kegg_barplot,
            "detected_columns": detected,
            "filtered_count": len(filtered_deg),
            "up_count": len(up_genes),
            "down_count": len(down_genes),
            "candidate_count": len(candidates),
        }
