"""Print basic statistics for the downloaded GenoPHI interaction matrix."""

from pathlib import Path

from phage_cocktail.analysis import summarize_interactions
from phage_cocktail.data import load_interactions

DATA_PATH = Path("data/raw/genophi_ecoli_subset.csv")


def main() -> None:
    if not DATA_PATH.exists():
        raise SystemExit(
            "Dataset not found. Run: python scripts/fetch_genophi.py"
        )

    frame = load_interactions(DATA_PATH)
    summary = summarize_interactions(frame)

    print("GenoPHI E. coli subset")
    print("----------------------")
    print(f"Interactions:      {summary['interactions']:,}")
    print(f"Bacterial strains: {summary['bacteria']:,}")
    print(f"Phages:            {summary['phages']:,}")
    print(f"Positive:          {summary['positive']:,}")
    print(f"Negative:          {summary['negative']:,}")
    print(f"Positive rate:     {summary['positive_rate']:.1%}")
    print(f"Observed matrix:   {summary['observed_fraction']:.1%}")


if __name__ == "__main__":
    main()
