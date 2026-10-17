"""Print basic statistics for the registered GenoPHI interaction matrix."""

from pathlib import Path

from phage_cocktail.registry import verify_registered_dataset

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "data/registry/genophi_ecoli_subset.v1.json"


def main() -> None:
    verified = verify_registered_dataset(MANIFEST_PATH)
    summary = verified.summary

    print(f"Dataset: {verified.manifest['dataset_id']}")
    print(f"SHA-256: {verified.sha256}")
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
