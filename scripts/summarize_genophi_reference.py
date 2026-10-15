"""Summarize metrics produced by the external GenoPHI reference run."""

from pathlib import Path

from phage_cocktail.genophi_reference import (
    find_performance_file,
    summarize_performance,
)

OUTPUT_DIR = Path("artifacts/genophi_demo")


def main() -> None:
    performance_file = find_performance_file(OUTPUT_DIR)
    summary = summarize_performance(performance_file)

    print(f"Source: {performance_file}")
    print(summary.round(4).to_string())


if __name__ == "__main__":
    main()
