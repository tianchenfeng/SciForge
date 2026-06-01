from __future__ import annotations

import argparse

from modules.module_silage_eval.service import SilageEvalService


def main():
    parser = argparse.ArgumentParser(description="Evaluate silage quality using entropy weight and membership score.")
    parser.add_argument("input_file")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    result = SilageEvalService.run(args.input_file, args.output_dir, logger=print)
    print(result["excel"])
    print(result["word"])


if __name__ == "__main__":
    main()
