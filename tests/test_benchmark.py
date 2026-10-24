import numpy as np
import pandas as pd
import pytest

from phage_cocktail.benchmark import classification_metrics, pairwise_baseline_predictions


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

    predictions = pairwise_baseline_predictions(train, test)

    assert classification_metrics(test["infection"], predictions["phage_prevalence"])["auroc"] == 1.0
    assert classification_metrics(test["infection"], predictions["coin_flip_0.5"])["auroc"] == 0.5
    # Beta(1,1) smoothing: (3 + 1) / (3 + 2) for the broad phage.
    assert np.allclose(predictions["phage_prevalence"], [0.8, 0.2, 0.8, 0.2])


def test_pairwise_baselines_reject_train_test_strain_overlap():
    train = _frame(["B1", "B2", "B3"])
    test = _frame(["B3", "B4"])

    with pytest.raises(ValueError, match="must be disjoint"):
        pairwise_baseline_predictions(train, test)
