"""Descriptive statistics for phage-host interaction datasets."""

from typing import Any

import pandas as pd

from .data import validate_interactions


def summarize_interactions(frame: pd.DataFrame) -> dict[str, Any]:
    """Summarize coverage, class balance, and matrix completeness."""
    validate_interactions(frame)

    n_interactions = len(frame)
    n_bacteria = frame["bacterium"].nunique()
    n_phages = frame["phage"].nunique()
    n_positive = int(frame["infection"].sum())
    n_negative = n_interactions - n_positive
    possible_pairs = n_bacteria * n_phages

    return {
        "interactions": n_interactions,
        "bacteria": n_bacteria,
        "phages": n_phages,
        "positive": n_positive,
        "negative": n_negative,
        "positive_rate": n_positive / n_interactions if n_interactions else 0.0,
        "possible_pairs": possible_pairs,
        "observed_fraction": (
            n_interactions / possible_pairs if possible_pairs else 0.0
        ),
    }
