"""Run leakage-safe baselines on the registered GenoPHI E. coli subset."""

from pathlib import Path

from phage_cocktail.benchmark import (
    evaluate_broadest_cocktails,
    evaluate_pairwise_baselines,
)
from phage_cocktail.data import load_interactions
from phage_cocktail.registry import verify_registered_dataset
from phage_cocktail.split import split_by_bacterium

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "data/registry/genophi_ecoli_subset.v1.json"


def main() -> None:
    verified = verify_registered_dataset(MANIFEST_PATH)
    frame = load_interactions(verified.path)
    assigned = split_by_bacterium(frame, seed=42)

    print(f"Dataset: {verified.manifest['dataset_id']}")
    print(f"SHA-256: {verified.sha256}")
    strain_counts = assigned.groupby("split")["bacterium"].nunique()
    interaction_counts = assigned["split"].value_counts()

    print("\nHeld-out-strain split")
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
