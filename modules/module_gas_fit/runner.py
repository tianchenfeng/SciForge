from __future__ import annotations

import argparse

from modules.module_gas_fit.service import GasFitService


def main():
    parser = argparse.ArgumentParser(description="Batch fit in vitro gas production curves.")
    parser.add_argument("input_file")
    parser.add_argument("output_dir")
    parser.add_argument("--model", choices=["exponential", "lag"], default="exponential")
    args = parser.parse_args()
    result = GasFitService.run(args.input_file, args.output_dir, args.model, logger=print)
    print(result["excel"])


if __name__ == "__main__":
    main()
