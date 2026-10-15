"""Naive benchmarks that future host-aware models must beat."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    log_loss,
    roc_auc_score,
)

from .data import validate_interactions


def _metrics(y_true: pd.Series, probability: np.ndarray) -> dict[str, float]:
    """Return ranking and probability-quality metrics."""
    y = y_true.to_numpy(dtype=int)
    p = np.asarray(probability, dtype=float)

    return {
        "auroc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "brier": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
    }


def evaluate_pairwise_baselines(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    smoothing: float = 1.0,
) -> pd.DataFrame:
    """Evaluate three host-agnostic probability baselines on held-out strains."""
    validate_interactions(train)
    validate_interactions(test)

    global_rate = float(train["infection"].mean())

    phage_stats = train.groupby("phage")["infection"].agg(["sum", "count"])
    phage_probability = (
        (phage_stats["sum"] + smoothing)
        / (phage_stats["count"] + 2.0 * smoothing)
    ).to_dict()

    predictions = {
        "coin_flip_0.5": np.full(len(test), 0.5),
        "global_prevalence": np.full(len(test), global_rate),
        "phage_prevalence": np.asarray(
            [phage_probability.get(phage, global_rate) for phage in test["phage"]],
            dtype=float,
        ),
    }

    rows = []
    for name, probability in predictions.items():
        rows.append({"baseline": name, **_metrics(test["infection"], probability)})

    return pd.DataFrame(rows).set_index("baseline")


def rank_phages_by_training_prevalence(train: pd.DataFrame) -> pd.DataFrame:
    """Rank candidate phages using only training-strain infection outcomes."""
    validate_interactions(train)

    ranking = (
        train.groupby("phage")["infection"]
        .agg(infected_train_strains="sum", train_strains="count", hit_rate="mean")
        .reset_index()
        .sort_values(
            ["hit_rate", "infected_train_strains", "phage"],
            ascending=[False, False, True],
        )
        .reset_index(drop=True)
    )
    return ranking


def evaluate_broadest_cocktails(
    train: pd.DataFrame,
    test: pd.DataFrame,
    cocktail_sizes: Iterable[int] = range(1, 6),
) -> pd.DataFrame:
    """Evaluate cocktails made by selecting the broadest training-set phages.

    The actual coverage is the fraction of held-out bacterial strains infected
    by at least one selected phage. The independence estimate is deliberately
    naive and is included to expose overlap/dependence between phages.
    """
    validate_interactions(train)
    validate_interactions(test)

    ranking = rank_phages_by_training_prevalence(train)
    rows: list[dict[str, object]] = []

    for size in cocktail_sizes:
        if size < 1 or size > len(ranking):
            raise ValueError("cocktail size must be between 1 and the number of phages")

        selected = ranking.head(size)
        names = selected["phage"].tolist()

        held_out = test[test["phage"].isin(names)]
        infected_by_strain = held_out.groupby("bacterium")["infection"].max()
        actual_coverage = float(infected_by_strain.mean())

        rates = selected["hit_rate"].to_numpy(dtype=float)
        independence_estimate = float(1.0 - np.prod(1.0 - rates))

        rows.append(
            {
                "cocktail_size": size,
                "phages": ", ".join(names),
                "train_independence_estimate": independence_estimate,
                "held_out_actual_coverage": actual_coverage,
            }
        )

    return pd.DataFrame(rows).set_index("cocktail_size")
