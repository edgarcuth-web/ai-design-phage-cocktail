import pandas as pd

from phage_cocktail.split import split_by_bacterium


def _complete_matrix(n_bacteria=20, n_phages=3):
    return pd.DataFrame(
        [
            {
                "bacterium": f"B{b:02d}",
                "phage": f"P{p:02d}",
                "infection": (b + p) % 2,
            }
            for b in range(n_bacteria)
            for p in range(n_phages)
        ]
    )


def test_split_keeps_each_bacterium_in_exactly_one_partition():
    result = split_by_bacterium(_complete_matrix(), seed=42)

    per_bacterium = result.groupby("bacterium")["split"].nunique()

    assert per_bacterium.max() == 1
    assert set(result["split"]) == {"train", "validation", "test"}


def test_split_is_stable_under_row_reordering():
    frame = _complete_matrix()

    a = split_by_bacterium(frame, seed=42)
    b = split_by_bacterium(frame.sample(frac=1, random_state=7), seed=42)

    map_a = a.drop_duplicates("bacterium").set_index("bacterium")["split"].to_dict()
    map_b = b.drop_duplicates("bacterium").set_index("bacterium")["split"].to_dict()

    assert map_a == map_b


def test_different_seed_changes_assignment():
    frame = _complete_matrix()

    a = split_by_bacterium(frame, seed=42)
    b = split_by_bacterium(frame, seed=43)

    map_a = a.drop_duplicates("bacterium").set_index("bacterium")["split"].to_dict()
    map_b = b.drop_duplicates("bacterium").set_index("bacterium")["split"].to_dict()

    assert map_a != map_b
