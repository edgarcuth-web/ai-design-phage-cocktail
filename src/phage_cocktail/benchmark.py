"""Host-agnostic comparator baselines and shared classification metrics."""

from __future__ import annotations

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
    """Score test rows with host-agnostic baselines fit on other strains' rows.

    Phage identities may occur in both partitions, so this measures unseen
    strains with known phages, not unseen phages.
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
