from __future__ import annotations

import argparse
import json
import sys

from modules.module_metabolomics.service import AnalysisConfig, MetabolomicsAnalysisService


def main(argv=None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="SciForge metabolomics runner")
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)
    with open(args.config, "r", encoding="utf-8") as handle: config = AnalysisConfig(**json.load(handle))
    def progress(stage, percent, message):
        print("PROGRESS\t" + json.dumps({"stage": stage, "percent": percent, "message": message}, ensure_ascii=False), flush=True)
    result = MetabolomicsAnalysisService.run(config, progress=progress)
    print("RESULT\t" + json.dumps(result, ensure_ascii=False), flush=True)
    return 0 if result.get("success") else 1


if __name__ == "__main__": sys.exit(main())
