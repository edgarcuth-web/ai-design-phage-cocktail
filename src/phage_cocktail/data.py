"""Utilities for loading and validating observed phage-host interaction rows.

Rows represent explicitly observed source labels. Missing bacterium-phage pairs
remain unknown; they are never filled in as negative outcomes by this module.
"""

from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ("bacterium", "phage", "infection")
COLUMN_ALIASES = {
    "strain": "bacterium",
    "interaction": "infection",
}


def normalize_interaction_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with common external column names mapped to our schema.

    GenoPHI uses strain, phage, interaction while this project internally uses
    bacterium, phage, infection. This is a column-name normalization only; it
    does not reinterpret the source phenotype or map it to a mechanistic target.
    """
    frame = frame.copy()
    rename: dict[str, str] = {}

    for source, target in COLUMN_ALIASES.items():
        if source in frame.columns:
            if target in frame.columns:
                raise ValueError(
                    f"ambiguous columns: both {source!r} and {target!r} are present"
                )
            rename[source] = target

    return frame.rename(columns=rename)


def validate_interactions(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate explicitly observed binary interaction rows without imputing."""
    missing = set(REQUIRED_COLUMNS).difference(frame.columns)
    if missing:
        names = ", ".join(sorted(missing))
        raise ValueError(f"missing required columns: {names}")

    required = frame.loc[:, list(REQUIRED_COLUMNS)]
    if required.isna().any().any():
        raise ValueError("required interaction fields cannot contain missing values")

    for column in ("bacterium", "phage"):
        if required[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"{column} identifiers cannot be empty")

    labels = set(required["infection"].unique())
    if not labels.issubset({0, 1}):
        raise ValueError("infection labels must be binary (0 or 1)")

    duplicated = required.duplicated(subset=["bacterium", "phage"], keep=False)
    if duplicated.any():
        raise ValueError("duplicate bacterium-phage pairs are not allowed")

    return frame


def load_interactions(path: str | Path) -> pd.DataFrame:
    """Load a CSV interaction table, normalize aliases, and validate it."""
    frame = pd.read_csv(path)
    frame = normalize_interaction_columns(frame)
    return validate_interactions(frame)
