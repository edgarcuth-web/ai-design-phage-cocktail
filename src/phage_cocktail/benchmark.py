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


def classification_metrics(y_true: pd.Series, probability: np.ndarray) -> dict[str, float]:
    """Return ranking and probability-quality metrics."""
    y = y_true.to_numpy(dtype=int)
    p = np.asarray(probability, dtype=float)

    return {
        "auroc": float(roc_auc_score(y, p)),
        "average_precision": float(average_precision_score(y, p)),
        "brier": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, p, labels=[0, 1])),
    }


def _require_disjoint_bacterial_strains(
    train: pd.DataFrame,
    test: pd.DataFrame,
) -> None:
    """Fail closed if a held-out-strain benchmark reuses target strains."""
    if train.empty or test.empty:
        raise ValueError("train and test partitions must both contain observations.")

    train_bacteria = set(train["bacterium"].astype(str))
    test_bacteria = set(test["bacterium"].astype(str))
    overlap = sorted(train_bacteria.intersection(test_bacteria))
    if overlap:
        preview = ", ".join(overlap[:5])
        suffix = " ..." if len(overlap) > 5 else ""
        raise ValueError(
            "train and test bacterial strains must be disjoint; "
            f"overlap: {preview}{suffix}"
        )


def pairwise_baseline_predictions(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    smoothing: float = 1.0,
) -> dict[str, np.ndarray]:
    """Evaluate host-agnostic baselines on unseen bacterial strains.

    Phage identities may occur in both partitions. This evaluates generalization
    to unseen strains for known phages, not cold-start generalization to unseen
    phages. Only explicitly observed test pairs are scored; missing pairs are
    never imputed as negative labels.
    """
    validate_interactions(train)
    validate_interactions(test)
    _require_disjoint_bacterial_strains(train, test)

    global_rate = float(train["infection"].mean())

    phage_stats = train.groupby("phage")["infection"].agg(["sum", "count"])
    phage_probability = (
        (phage_stats["sum"] + smoothing)
        / (phage_stats["count"] + 2.0 * smoothing)
    ).to_dict()

    return {
        "coin_flip_0.5": np.full(len(test), 0.5),
        "global_prevalence": np.full(len(test), global_rate),
        "phage_prevalence": np.asarray(
            [phage_probability.get(phage, global_rate) for phage in test["phage"]],
            dtype=float,
        ),
    }


def evaluate_pairwise_baselines(
    train: pd.DataFrame, test: pd.DataFrame, *, smoothing: float = 1.0,
) -> pd.DataFrame:
    """Evaluate the shared baseline predictions on the same observed rows."""
    predictions = pairwise_baseline_predictions(train, test, smoothing=smoothing)
    rows = []
    for name, probability in predictions.items():
        rows.append({"baseline": name, **classification_metrics(test["infection"], probability)})

    return pd.DataFrame(rows).set_index("baseline")


def rank_phages_by_training_prevalence(train: pd.DataFrame) -> pd.DataFrame:
    """Rank candidate phages using only explicitly observed training outcomes."""
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


def _require_complete_cocktail_matrix(
    test: pd.DataFrame,
    selected_phages: list[str],
) -> None:
    """Require observations for every selected-phage/held-out-strain pair."""
    bacteria = sorted(test["bacterium"].astype(str).unique())
    observed = set(
        zip(test["bacterium"].astype(str), test["phage"].astype(str), strict=True)
    )
    missing = [
        (bacterium, phage)
        for bacterium in bacteria
        for phage in selected_phages
        if (bacterium, phage) not in observed
    ]
    if missing:
        preview = ", ".join(f"({b}, {p})" for b, p in missing[:5])
        suffix = " ..." if len(missing) > 5 else ""
        raise ValueError(
            "cannot calculate whole-cohort cocktail coverage: selected-phage "
            "pairs are missing from the test matrix; missing pairs are unknown, "
            f"not negative. Missing: {preview}{suffix}"
        )


def evaluate_broadest_cocktails(
    train: pd.DataFrame,
    test: pd.DataFrame,
    cocktail_sizes: Iterable[int] = range(1, 6),
) -> pd.DataFrame:
    """Evaluate offline union coverage selected by training-set phage prevalence.

    This is a union of single-phage labels, not a measured mixture outcome.
    Held-out coverage is only reported when every selected phage has an
    observed outcome for every test strain. Missing pairs are unknown, not
    negative, and cannot be silently omitted from the coverage denominator.
    The independence estimate is deliberately naive and exposes overlap or
    dependence between phages.
    """
    validate_interactions(train)
    validate_interactions(test)
    _require_disjoint_bacterial_strains(train, test)

    ranking = rank_phages_by_training_prevalence(train)
    rows: list[dict[str, object]] = []

    for size in cocktail_sizes:
        if size < 1 or size > len(ranking):
            raise ValueError("cocktail size must be between 1 and the number of phages")

        selected = ranking.head(size)
        names = selected["phage"].astype(str).tolist()
        _require_complete_cocktail_matrix(test, names)

        held_out = test[test["phage"].astype(str).isin(names)]
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
