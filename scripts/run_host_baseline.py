"""Audit locked real data, freeze genomic groups, and run the fixed host baseline."""
import argparse
from pathlib import Path

from phage_cocktail.host_study import run_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True, help="A new run directory")
    parser.add_argument("--fetch", action="store_true", help="Download missing locked public resources")
    args = parser.parse_args()
    run_study(Path(__file__).resolve().parents[1], args.data_root, args.outdir, fetch=args.fetch)


if __name__ == "__main__":
    main()
