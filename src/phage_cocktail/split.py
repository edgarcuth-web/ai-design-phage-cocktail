"""Leakage-safe train/validation/test splitting by bacterial strain."""

from __future__ import annotations

import pandas as pd

from .data import validate_interactions


def _fnv1a_32(text: str) -> int:
    """Return a stable 32-bit FNV-1a hash for deterministic ordering."""
    value = 2166136261
    for byte in text.encode("utf-8"):
        value ^= byte
        value = (value * 16777619) & 0xFFFFFFFF
    return value


def _split_counts(n: int, fractions: tuple[float, float, float]) -> tuple[int, int, int]:
    """Allocate n groups using the largest-remainder method."""
    if n < 3:
        raise ValueError("at least three bacterial strains are required")

    if any(f <= 0 for f in fractions):
        raise ValueError("all split fractions must be positive")
    if abs(sum(fractions) - 1.0) > 1e-9:
        raise ValueError("split fractions must sum to 1")

    raw = [n * f for f in fractions]
    counts = [int(x) for x in raw]
    remaining = n - sum(counts)

    order = sorted(
        range(3),
        key=lambda i: (raw[i] - counts[i], -i),
        reverse=True,
    )
    for i in order[:remaining]:
        counts[i] += 1

    return counts[0], counts[1], counts[2]


def split_by_bacterium(
    frame: pd.DataFrame,
    *,
    train_fraction: float = 0.70,
    validation_fraction: float = 0.15,
    test_fraction: float = 0.15,
    seed: int = 42,
) -> pd.DataFrame:
    """Assign every bacterial strain wholly to train, validation, or test.

    The assignment is stable across row order and pandas/sklearn versions.
    No interaction from a held-out strain can appear in training.
    """
    validate_interactions(frame)

    bacteria = sorted(
        (str(x) for x in frame["bacterium"].unique()),
        key=lambda name: (_fnv1a_32(f"{seed}:{name}"), name),
    )
    n_train, n_validation, _ = _split_counts(
        len(bacteria),
        (train_fraction, validation_fraction, test_fraction),
    )

    assignments: dict[str, str] = {}
    for name in bacteria[:n_train]:
        assignments[name] = "train"
    for name in bacteria[n_train : n_train + n_validation]:
        assignments[name] = "validation"
    for name in bacteria[n_train + n_validation :]:
        assignments[name] = "test"

    result = frame.copy()
    result["bacterium"] = result["bacterium"].astype(str)
    result["split"] = result["bacterium"].map(assignments)

    if result["split"].isna().any():
        raise RuntimeError("failed to assign one or more bacterial strains")

    return result


def split_by_host_group(frame: pd.DataFrame, groups: dict[str, str], *, seed: int = 42,
                        fractions=(0.70, 0.15, 0.15)) -> pd.DataFrame:
    """Assign whole host groups, preserving the legacy split implementation."""
    validate_interactions(frame)
    hosts = set(frame.bacterium.astype(str))
    if hosts - set(groups) or any(not isinstance(groups[h], str) or not groups[h] for h in hosts):
        raise ValueError("every bacterial strain requires a non-empty host group")
    representatives = sorted(set(groups[h] for h in hosts),
                             key=lambda group: (_fnv1a_32(f"{seed}:{group}"), group))
    n_train, n_validation, _ = _split_counts(len(representatives), tuple(fractions))
    assignment = {g: ("train" if i < n_train else "validation" if i < n_train + n_validation else "test")
                  for i, g in enumerate(representatives)}
    result = frame.copy()
    result["host_group"] = result.bacterium.map(groups)
    result["split"] = result.host_group.map(assignment)
    return result
