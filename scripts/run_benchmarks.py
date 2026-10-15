"""Run leakage-safe baselines on the pinned GenoPHI E. coli subset."""

from pathlib import Path

from phage_cocktail.benchmark import (
    evaluate_broadest_cocktails,
    evaluate_pairwise_baselines,
)
from phage_cocktail.data import load_interactions
from phage_cocktail.split import split_by_bacterium

DATA_PATH = Path("data/raw/genophi_ecoli_subset.csv")


def main() -> None:
    if not DATA_PATH.exists():
        raise SystemExit("Dataset not found. Run: python scripts/fetch_genophi.py")

    frame = load_interactions(DATA_PATH)
    assigned = split_by_bacterium(frame, seed=42)

    strain_counts = assigned.groupby("split")["bacterium"].nunique()
    interaction_counts = assigned["split"].value_counts()

    print("Held-out-strain split")
    print("---------------------")
    for split in ("train", "validation", "test"):
        print(
            f"{split:10s}: "
            f"{strain_counts[split]:3d} strains, "
            f"{interaction_counts[split]:5d} interactions"
        )

    train = assigned[assigned["split"] == "train"].drop(columns="split")
    validation = assigned[assigned["split"] == "validation"].drop(columns="split")
    test = assigned[assigned["split"] == "test"].drop(columns="split")

    # Validation is reserved for later model selection. The naive benchmark is
    # reported once on test without tuning against it.
    print("\nPairwise baselines on held-out TEST strains")
    print("------------------------------------------")
    print(evaluate_pairwise_baselines(train, test).round(4).to_string())

    print("\nBroadest-phage cocktails on held-out TEST strains")
    print("------------------------------------------------")
    print(evaluate_broadest_cocktails(train, test).round(4).to_string())

    print(
        "\nValidation interactions reserved for later model selection:",
        f"{len(validation):,}",
    )


if __name__ == "__main__":
    main()
