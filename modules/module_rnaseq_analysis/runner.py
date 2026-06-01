from __future__ import annotations

import argparse

from modules.module_rnaseq_analysis.service import RNASeqAnalysisService


def main():
    parser = argparse.ArgumentParser(description="Comprehensive RNASeq result analysis.")
    parser.add_argument("deg_file")
    parser.add_argument("output_dir")
    parser.add_argument("--log2fc", type=float, default=1.0)
    parser.add_argument("--padj", type=float, default=0.05)
    parser.add_argument("--keywords", default="")
    parser.add_argument("--annotation")
    parser.add_argument("--go")
    parser.add_argument("--kegg")
    parser.add_argument("--expression")
    args = parser.parse_args()
    result = RNASeqAnalysisService.run(
        args.deg_file,
        args.output_dir,
        args.log2fc,
        args.padj,
        args.keywords,
        args.annotation,
        args.go,
        args.kegg,
        args.expression,
        logger=print,
    )
    print(result["excel"])
    print(result["word"])


if __name__ == "__main__":
    main()
