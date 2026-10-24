"""Small additive pair baselines and paired group-level evaluation."""
from dataclasses import dataclass
from itertools import product
import warnings

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss, roc_auc_score

from .host_baseline import HostFeatureEncoder

COMPLEMENT = str.maketrans("ACGT", "TGCA")
FOURMERS = sorted({min(s, s.translate(COMPLEMENT)[::-1])
                   for s in map("".join, product("ACGT", repeat=4))})
FOURMER_INDEX = {s: i for i, s in enumerate(FOURMERS)}


def canonical_fourmers(payload):
    """Count admitted DNA windows within contigs; ambiguous windows are skipped."""
    contigs, pieces = [], None
    for line in payload.decode("ascii").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if pieces is not None:
                contigs.append("".join(pieces).upper())
            pieces = []
        elif pieces is None:
            raise ValueError("DNA FASTA requires a header")
        else:
            pieces.append(line)
    if pieces is not None:
        contigs.append("".join(pieces).upper())
    if not contigs or any(not s or set(s) - set("ACGTRYSWKMBDHVN") for s in contigs):
        raise ValueError("invalid or empty DNA contig")
    result = np.zeros(len(FOURMERS), dtype=float)
    for sequence in contigs:
        for start in range(len(sequence) - 3):
            word = sequence[start:start + 4]
            if set(word) <= set("ACGT"):
                result[FOURMER_INDEX[min(word, word.translate(COMPLEMENT)[::-1])]] += 1
    if not result.sum():
        raise ValueError("sequence has no valid 4-mer windows")
    return result / result.sum()


@dataclass
class PairEncoder:
    host: HostFeatureEncoder
    phage_columns: list
    phage_means: list
    phage_scales: list

    @classmethod
    def fit(cls, train, hosts, phages, numeric, categorical, phage_digests):
        host_ids, phage_ids = sorted(set(train.bacterium)), sorted(set(train.phage))
        encoder = HostFeatureEncoder.fit(hosts.loc[host_ids], numeric, categorical)
        unique = {}
        for name in phage_ids:
            digest = phage_digests[name]
            if digest in unique and not phages.loc[name].equals(phages.loc[unique[digest]]):
                raise ValueError("identical phage assemblies have inconsistent features")
            unique.setdefault(digest, name)
        values = phages.loc[list(unique.values())].to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("phage features must be finite")
        scales = values.std(axis=0)
        scales[scales == 0] = 1
        return cls(encoder, phages.columns.tolist(), values.mean(axis=0).tolist(), scales.tolist())

    def transform(self, pairs, hosts, phages, kind):
        blocks = []
        if kind in {"host", "pair"}:
            blocks.append(self.host.transform(hosts.loc[pairs.bacterium].reset_index(drop=True)))
        if kind in {"phage", "pair"}:
            values = phages.loc[pairs.phage, self.phage_columns].to_numpy(dtype=float)
            blocks.append((values - self.phage_means) / self.phage_scales)
        if not blocks or not all(np.isfinite(b).all() for b in blocks):
            raise ValueError("unknown representation or unsupported/non-finite feature input")
        return np.column_stack(blocks)

    def to_dict(self):
        return {"host": self.host.to_dict(), "phage_columns": self.phage_columns,
                "phage_means": self.phage_means, "phage_scales": self.phage_scales,
                "statistics_scope": "unique fit hosts and unique fit phage assembly digests"}

    @classmethod
    def from_dict(cls, data):
        return cls(HostFeatureEncoder(**data["host"]), data["phage_columns"],
                   data["phage_means"], data["phage_scales"])


def fit_logistic(x, y, weights, *, C, max_iter):
    if not np.isfinite(x).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("fit requires finite features and positive weights")
    if set(np.unique(y)) != {0, 1}:
        raise ValueError("fit requires both source-label classes")
    model = LogisticRegression(C=C, solver="lbfgs", max_iter=max_iter, tol=0.0001,
                               fit_intercept=True, class_weight=None)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        model.fit(x, y, sample_weight=weights)
    return {"C": C, "coefficients": model.coef_[0].tolist(), "intercept": float(model.intercept_[0]),
            "iterations": int(model.n_iter_[0]), "score_semantics": "uncalibrated_source_mlc_positive"}


def predict_logistic(model, x):
    # Batch BLAS kernels can split exact feature ties by a few ULPs. AP depends
    # on those ties, so accumulate each row in the same feature order.
    x, coefficients = np.asarray(x, dtype=float), np.asarray(model["coefficients"])
    if x.ndim != 2 or x.shape[1] != len(coefficients) or not np.isfinite(x).all():
        raise ValueError("prediction matrix must match finite model inputs")
    logits = np.zeros(len(x), dtype=float)
    for column, coefficient in enumerate(coefficients):
        logits += x[:, column] * coefficient
    return expit(logits + model["intercept"])


def sequence_pair_weights(rows, host_digests, phage_digests):
    keys = pd.MultiIndex.from_arrays([rows.bacterium.map(host_digests), rows.phage.map(phage_digests)])
    if keys.to_frame(index=False).isna().any().any():
        raise ValueError("every row needs admitted host and phage assembly identities")
    counts = pd.Series(list(keys)).value_counts()
    return np.asarray([1. / counts[key] for key in keys])


def _ap_function(y, scores):
    order = np.argsort(-np.asarray(scores), kind="stable")
    labels, sorted_scores = np.asarray(y)[order], np.asarray(scores)[order]
    ends = np.r_[np.flatnonzero(np.diff(sorted_scores)), len(order) - 1]

    def calculate(weights):
        w = np.asarray(weights, dtype=float)[order]
        tp, totals = np.cumsum(w * labels)[ends], np.cumsum(w)[ends]
        if tp[-1] <= 0:
            raise ValueError("AP requires weighted positive observations")
        precision = np.divide(tp, totals, out=np.zeros_like(tp), where=totals > 0)
        return float(np.sum(np.diff(np.r_[0., tp]) * precision) / tp[-1])
    return calculate


def weighted_ap(y, scores, weights):
    return _ap_function(y, scores)(weights)


def metric_report(rows, scores, weights):
    """Report support explicitly; undefined class/group metrics remain null."""
    y, scores, weights = rows.infection.to_numpy(), np.asarray(scores), np.asarray(weights)
    keep = np.isfinite(scores) & (weights > 0)
    result = {"rows": len(rows), "scored_rows": int(keep.sum()), "coverage": float(keep.mean()),
              "positive_rows": int(y.sum())}
    if not keep.any():
        return {**result, "status": "not_estimated_unsupported_phages"}
    y, p, w = y[keep], scores[keep], weights[keep]
    both = len(np.unique(y)) == 2
    result.update(status="uncalibrated_score", average_precision=weighted_ap(y, p, w) if both else None,
                  auroc=float(roc_auc_score(y, p, sample_weight=w)) if both else None,
                  brier=float(np.average((p - y) ** 2, weights=w)),
                  log_loss=float(log_loss(y, p, labels=[0, 1], sample_weight=w)),
                  weighted_total=float(w.sum()), weighted_positive=float(w @ y))
    admitted = rows.loc[keep].reset_index(drop=True)
    for column in ("bacterium", "phage"):
        values, undefined = [], 0
        for positions in admitted.groupby(column, sort=True).indices.values():
            if len(np.unique(y[positions])) != 2:
                undefined += 1
            else:
                values.append(weighted_ap(y[positions], p[positions], w[positions]))
        result[f"macro_{column}_ap"] = float(np.mean(values)) if values else None
        result[f"macro_{column}_defined"] = len(values)
        result[f"macro_{column}_undefined"] = undefined
    return result


def grouped_ap_contrast(rows, scores, comparator, weights, regime, *, resamples, seed, minimum_groups=5):
    """Resample whole groups, multiplying multiplicities for crossed axes."""
    axes = {"cold_host": ["host_group"], "cold_phage": ["phage_group"],
            "double_cold": ["host_group", "phage_group"], "warm_pair": ["host_group", "phage_group"]}[regime]
    y, weights = rows.infection.to_numpy(), np.asarray(weights)
    evaluate, baseline = _ap_function(y, scores), _ap_function(y, comparator)
    result = {"delta_average_precision": evaluate(weights) - baseline(weights),
              "interval": None, "improvement_supported": False, "resamples": resamples, "seed": seed}
    encoded = [np.unique(rows[axis].to_numpy(), return_inverse=True) for axis in axes]
    result["resampled_groups"] = {axis: len(groups) for axis, (groups, _) in zip(axes, encoded)}
    if any(len(groups) < minimum_groups for groups, _ in encoded):
        return {**result, "status": "descriptive_insufficient_groups", "valid_resamples": 0}
    rng, deltas, skipped = np.random.default_rng(seed), [], 0
    for _ in range(resamples):
        current = weights.copy()
        for groups, indices in encoded:
            count = len(groups)
            multiplicity = np.bincount(rng.integers(count, size=count), minlength=count)
            current *= multiplicity[indices]
        positive, total = current @ y, current.sum()
        if positive <= 0 or positive >= total:
            skipped += 1
            continue
        deltas.append(evaluate(current) - baseline(current))
    result.update(valid_resamples=len(deltas), single_class_resamples=skipped)
    if not deltas:
        return {**result, "status": "descriptive_no_valid_resamples"}
    lower, upper = np.quantile(deltas, [0.025, 0.975]).tolist()
    return {**result, "status": "conditional_development_interval", "interval": [lower, upper],
            "improvement_supported": lower > 0}
