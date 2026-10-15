"""Simple baselines for combining phage infection probabilities."""

from collections.abc import Iterable

import numpy as np


def cocktail_success_probability(probabilities: Iterable[float]) -> float:
    """Estimate P(at least one phage infects) under independence.

    Parameters
    ----------
    probabilities:
        Per-phage infection probabilities, each in the closed interval [0, 1].

    Returns
    -------
    float
        Probability that at least one phage infects.

    Notes
    -----
    This is a deliberately naive baseline. Real phages can be dependent because
    they may share receptors, defence-system sensitivities, or interact with one
    another.
    """
    values = np.asarray(list(probabilities), dtype=float)

    if values.ndim != 1:
        raise ValueError("probabilities must be a one-dimensional sequence")
    if values.size == 0:
        raise ValueError("probabilities must contain at least one value")
    if not np.all(np.isfinite(values)):
        raise ValueError("probabilities must all be finite")
    if np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("probabilities must lie between 0 and 1")

    return float(1.0 - np.prod(1.0 - values))
