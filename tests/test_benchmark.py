import pandas as pd

from phage_cocktail.benchmark import (
    evaluate_broadest_cocktails,
    evaluate_pairwise_baselines,
)


def _frame(bacteria):
    rows = []
    for bacterium in bacteria:
        rows.extend(
            [
                {"bacterium": bacterium, "phage": "broad", "infection": 1},
                {"bacterium": bacterium, "phage": "narrow", "infection": 0},
            ]
        )
    return pd.DataFrame(rows)


def test_phage_prevalence_discriminates_broad_from_narrow_phage():
    train = _frame(["B1", "B2", "B3"])
    test = _frame(["B4", "B5"])

    result = evaluate_pairwise_baselines(train, test)

    assert result.loc["phage_prevalence", "auroc"] == 1.0
    assert result.loc["coin_flip_0.5", "auroc"] == 0.5


def test_broadest_cocktail_selects_broad_phage_first():
    train = _frame(["B1", "B2", "B3"])
    test = _frame(["B4", "B5"])

    result = evaluate_broadest_cocktails(train, test, cocktail_sizes=[1])

    assert result.loc[1, "phages"] == "broad"
    assert result.loc[1, "held_out_actual_coverage"] == 1.0
