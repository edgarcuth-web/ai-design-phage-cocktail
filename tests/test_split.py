import pandas as pd

from phage_cocktail.split import assign_partitions, split_by_host_group

NAMES = [f"G{i:02d}" for i in range(20)]


def test_every_partition_is_used_with_largest_remainder_counts():
    assignment = assign_partitions(NAMES, seed=42)

    counts = pd.Series(assignment).value_counts()
    assert counts.to_dict() == {"train": 14, "validation": 3, "test": 3}


def test_assignment_is_stable_under_input_order_and_changes_with_seed():
    assert assign_partitions(NAMES, seed=42) == assign_partitions(reversed(NAMES), seed=42)
    assert assign_partitions(NAMES, seed=42) != assign_partitions(NAMES, seed=43)


def test_host_groups_are_never_split_across_partitions():
    frame = pd.DataFrame(
        [
            {"bacterium": f"B{b:02d}", "phage": f"P{p}", "infection": (b + p) % 2}
            for b in range(40)
            for p in range(3)
        ]
    )
    groups = {f"B{b:02d}": f"G{b // 2:02d}" for b in range(40)}

    result = split_by_host_group(frame, groups, seed=42)

    assert result.groupby("host_group")["split"].nunique().max() == 1
    assert set(result["split"]) == {"train", "validation", "test"}
