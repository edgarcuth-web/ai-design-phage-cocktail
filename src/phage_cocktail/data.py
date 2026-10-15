"""Utilities for loading pairwise phage-host interaction tables."""

from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = {"bacterium", "phage", "infection"}


def load_interactions(path: str | Path) -> pd.DataFrame:
    """Load and validate a CSV interaction table."""
    frame = pd.read_csv(path)

    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        names = ", ".join(sorted(missing))
        raise ValueError(f"missing required columns: {names}")

    if frame[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValueError("required interaction fields cannot contain missing values")

    labels = set(frame["infection"].unique())
    if not labels.issubset({0, 1}):
        raise ValueError("infection labels must be binary (0 or 1)")

    return frame
