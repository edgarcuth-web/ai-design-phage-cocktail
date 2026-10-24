"""Small known-phage baseline with train-only host preprocessing and JSON weights."""
from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score

from .data import validate_interactions

MISSING = "__NOT_REPORTED__"


def _numeric(frame, columns):
    values = frame.loc[:, columns].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("defense counts must be finite and non-negative")
    return values


def _category(series):
    if series.dropna().astype(str).eq(MISSING).any():
        raise ValueError("source category collides with missingness sentinel")
    return series.fillna(MISSING).astype(str)


@dataclass
class HostFeatureEncoder:
    numeric_columns: list[str]
    categorical_columns: list[str]
    means: list[float]
    scales: list[float]
    categories: dict[str, list[str]]

    @classmethod
    def fit(cls, training_hosts, numeric_columns, categorical_columns):
        if training_hosts.empty or training_hosts.index.has_duplicates:
            raise ValueError("preprocessing requires unique, non-empty training hosts")
        numeric_columns, categorical_columns = list(numeric_columns), list(categorical_columns)
        values = _numeric(training_hosts, numeric_columns)
        means, scales = values.mean(axis=0), values.std(axis=0)
        scales[scales == 0] = 1
        categories = {
            column: sorted(_category(training_hosts[column]).unique().tolist())
            for column in categorical_columns
        }
        return cls(numeric_columns, categorical_columns, means.tolist(), scales.tolist(), categories)

    def transform(self, hosts):
        if hosts.index.has_duplicates:
            raise ValueError("duplicate host feature rows")
        pieces = [(_numeric(hosts, self.numeric_columns) - self.means) / self.scales]
        for column in self.categorical_columns:
            values = _category(hosts[column])
            indicators = [(values == value).to_numpy(dtype=float) for value in self.categories[column]]
            pieces.append(np.column_stack(indicators))
        return np.column_stack(pieces)

    def to_dict(self):
        return {
            "numeric_columns": self.numeric_columns,
            "categorical_columns": self.categorical_columns,
            "means": self.means,
            "scales": self.scales,
            "categories": self.categories,
        }


def fit_host_model(
    train, features, numeric_columns, categorical_columns, *, C, max_iter, sample_weight=None
):
    """Fit one fixed logistic model per training phage, using no held-out hosts."""
    validate_interactions(train)
    train = train.reset_index(drop=True)
    if sample_weight is None:
        weights = np.ones(len(train))
    else:
        weights = np.asarray(sample_weight, dtype=float)
    if weights.shape != (len(train),) or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError("sample weights must be aligned, finite and positive")
    train_hosts = sorted(train.bacterium.unique())
    if not train_hosts or set(train_hosts) - set(features.index):
        raise ValueError("missing training host features")
    encoder = HostFeatureEncoder.fit(
        features.loc[train_hosts], numeric_columns, categorical_columns
    )
    host_matrix = encoder.transform(features.loc[train_hosts])
    row_index = {host: i for i, host in enumerate(train_hosts)}
    models = {}
    for phage, observations in train.groupby("phage", sort=True):
        y = observations.infection.to_numpy(dtype=int)
        local_weights = weights[observations.index]
        if len(np.unique(y)) < 2:
            models[phage] = {
                "kind": "single_class_prevalence",
                "score": float((local_weights @ y + 1) / (local_weights.sum() + 2)),
                "training_rows": len(y),
            }
            continue
        x = host_matrix[[row_index[host] for host in observations.bacterium]]
        classifier = LogisticRegression(C=C, solver="lbfgs", max_iter=max_iter)
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            classifier.fit(x, y, sample_weight=None if sample_weight is None else local_weights)
        models[phage] = {
            "kind": "logistic",
            "coefficients": classifier.coef_[0].tolist(),
            "intercept": float(classifier.intercept_[0]),
            "training_rows": len(y),
            "iterations": int(classifier.n_iter_[0]),
        }
    return {
        "schema_version": 1,
        "score_semantics": "uncalibrated_source_mlc_positive",
        "training_hosts": train_hosts,
        "preprocessing": encoder.to_dict(),
        "phages": models,
    }


def predict_host_model(model, pairs, features):
    """Predict from host features and known phage identity, never from labels."""
    if set(pairs.phage) - set(model["phages"]):
        raise ValueError("unsupported phage identity; this baseline has no cold-start model")
    hosts = sorted(pairs.bacterium.unique())
    if set(hosts) - set(features.index):
        raise ValueError("missing prediction host features")
    encoder = HostFeatureEncoder(**model["preprocessing"])
    matrix = encoder.transform(features.loc[hosts])
    indices = {host: i for i, host in enumerate(hosts)}
    scores = []
    for row in pairs.itertuples(index=False):
        estimator = model["phages"][row.phage]
        if estimator["kind"] == "single_class_prevalence":
            scores.append(estimator["score"])
        else:
            row_features = matrix[indices[row.bacterium]]
            z = float(row_features @ estimator["coefficients"] + estimator["intercept"])
            # Stable sigmoid, including extreme inputs outside the training range.
            scores.append(1 / (1 + np.exp(-z)) if z >= 0 else np.exp(z) / (1 + np.exp(z)))
    return np.asarray(scores, dtype=float)


def paired_group_bootstrap(
    pairs, host_scores, baseline_scores, *, resamples, seed, confidence=0.95
):
    """Paired percentile interval for AP gain, resampling complete host groups."""
    y = pairs.infection.to_numpy(dtype=int)
    a, b = np.asarray(host_scores), np.asarray(baseline_scores)
    if len(y) != len(a) or a.shape != b.shape or len(np.unique(y)) != 2:
        raise ValueError("bootstrap needs paired scores and both observed classes")
    if resamples < 2 or not 0 < confidence < 1 or pairs.host_group.isna().any():
        raise ValueError("invalid grouped bootstrap configuration")
    groups = sorted(pairs.host_group.unique())
    if len(groups) < 2:
        raise ValueError("bootstrap requires at least two held-out groups")
    indices = [np.flatnonzero(pairs.host_group.to_numpy() == group) for group in groups]
    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(resamples):
        drawn = rng.integers(len(groups), size=len(groups))
        sample = np.concatenate([indices[i] for i in drawn])
        if len(np.unique(y[sample])) < 2:
            continue
        deltas.append(
            float(
                average_precision_score(y[sample], a[sample])
                - average_precision_score(y[sample], b[sample])
            )
        )
    if not deltas:
        raise ValueError("all grouped bootstrap resamples are single-class")
    alpha = (1 - confidence) / 2
    interval = np.quantile(deltas, [alpha, 1 - alpha]).tolist()
    delta = float(average_precision_score(y, a) - average_precision_score(y, b))
    return {
        "delta_average_precision": delta,
        "interval": interval,
        "confidence_level": confidence,
        "groups": len(groups),
        "requested_resamples": resamples,
        "valid_resamples": len(deltas),
        "single_class_resamples": resamples - len(deltas),
        "seed": seed,
        "improvement_supported": bool(interval[0] > 0),
    }
