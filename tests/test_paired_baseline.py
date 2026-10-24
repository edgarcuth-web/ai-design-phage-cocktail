"""Scientific invariants for paired-model fitting and crossed evaluation."""
import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score

from phage_cocktail.paired_baseline import (
    canonical_fourmers, PairEncoder, fit_logistic, predict_logistic,
    sequence_pair_weights, grouped_ap_contrast, weighted_ap, metric_report,
)
from phage_cocktail.host_baseline import fit_host_model


def test_fourmers_are_strand_invariant_skip_ambiguity_and_never_bridge_contigs():
    a = canonical_fourmers(b">a\nAAAACG\n>b\nTTTT\n")
    b = canonical_fourmers(b">different\nCGTTTT\n>other\nAAAA\n")
    assert len(a) == 136 and np.isclose(a.sum(), 1)
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(canonical_fourmers(b">a\nAAAANNNNTTTT\n"),
                                  canonical_fourmers(b">a\nAAAA\n>b\nTTTT\n"))
    with pytest.raises(ValueError, match="valid 4-mer"):
        canonical_fourmers(b">a\nAAA\n>b\nA\n")
    with pytest.raises(ValueError, match="DNA|header"):
        canonical_fourmers(b"AAAA")


def fixtures():
    hosts = pd.DataFrame({"d": [0., 2., 1000.], "surface": ["A", "B", "HELD_ONLY"]}, index=["H1", "H2", "HX"])
    phages = pd.DataFrame({"f1": [0., 2., 1000.], "f2": [2., 0., 1000.]}, index=["P1", "P2", "PX"])
    train = pd.DataFrame({"bacterium": ["H1", "H1", "H2", "H2"],
                          "phage": ["P1", "P2", "P1", "P2"], "infection": [0, 1, 1, 0]})
    return hosts, phages, train


def test_encoders_fit_unique_training_entities_and_predictions_ignore_labels():
    hosts, phages, train = fixtures()
    encoder = PairEncoder.fit(train, hosts, phages, ["d"], ["surface"], {p: p for p in phages.index})
    assert encoder.host.means == [1.]
    assert encoder.phage_means == [1., 1.]
    assert "HELD_ONLY" not in encoder.host.categories["surface"]
    x = encoder.transform(train, hosts, phages, "pair")
    model = fit_logistic(x, train.infection.to_numpy(), np.ones(4), C=0.1, max_iter=2000)
    restored = PairEncoder.from_dict(encoder.to_dict())
    np.testing.assert_array_equal(x, restored.transform(train, hosts, phages, "pair"))
    np.testing.assert_array_equal(predict_logistic(model, x), predict_logistic(model, encoder.transform(
        train.assign(infection=1), hosts, phages, "pair")))
    held = train.iloc[:1].assign(bacterium="HX", phage="PX")
    assert np.isfinite(encoder.transform(held, hosts, phages, "pair")).all()


def test_sequence_identity_deduplicates_phage_statistics_and_row_weights():
    hosts, phages, train = fixtures()
    phages.loc["COPY"] = phages.loc["P1"]
    extra = train[train.phage == "P1"].assign(phage="COPY")
    rows = pd.concat([train, extra], ignore_index=True)
    digests = {"P1": "same", "COPY": "same", "P2": "other", "PX": "held"}
    encoder = PairEncoder.fit(rows, hosts, phages, ["d"], ["surface"], digests)
    assert encoder.phage_means == [1., 1.]
    weights = sequence_pair_weights(rows, {h: h for h in hosts.index}, digests)
    assert np.isclose(weights.sum(), 4)
    assert set(weights) == {0.5, 1.0}


def test_weighted_ap_matches_sklearn_with_ties_and_zero_weight_rows():
    y = np.array([0, 1, 1, 0, 1])
    scores = np.array([0.8, 0.8, 0.4, 0.4, 0.1])
    weights = np.array([0., 2., 1., 3., 0.5])
    assert np.isclose(weighted_ap(y, scores, weights), average_precision_score(y, scores, sample_weight=weights))
    integer_weights = np.array([1, 2, 0, 3, 1])
    assert np.isclose(weighted_ap(y, scores, integer_weights),
                      average_precision_score(y, scores, sample_weight=integer_weights))


def test_identical_features_keep_exact_ties_independent_of_batch_shape():
    from threadpoolctl import threadpool_limits
    rng = np.random.default_rng(1)
    base = rng.normal(size=(7, 257))
    indices = np.arange(503) % len(base)
    model = {"coefficients": rng.normal(scale=.1, size=257).tolist(), "intercept": .1}
    with threadpool_limits(limits=1):
        expected = predict_logistic(model, base)[indices]
        actual = predict_logistic(model, base[indices])
    np.testing.assert_array_equal(actual, expected)


def test_legacy_weighted_fallback_and_unsupported_coverage_are_explicit():
    hosts, _, train = fixtures()
    rows = train.assign(infection=0)
    model = fit_host_model(rows, hosts, ["d"], ["surface"], C=0.1, max_iter=2000,
                           sample_weight=np.full(4, 0.5))
    assert np.isclose(model["phages"]["P1"]["score"], 1/3)
    report = metric_report(train, np.full(4, np.nan), np.ones(4))
    assert report["coverage"] == 0
    assert report["status"] == "not_estimated_unsupported_phages"
    with pytest.raises(ValueError, match="weights"):
        fit_host_model(rows, hosts, ["d"], ["surface"], C=0.1, max_iter=2000, sample_weight=[1, 2])


def test_crossed_bootstrap_equal_predictions_have_zero_delta_and_report_small_panels():
    rows = pd.DataFrame({"host_group": np.repeat(list("ABCDE"), 5),
                         "phage_group": list("VWXYZ") * 5,
                         "infection": [0, 1, 0, 1, 1] * 5})
    scores = np.linspace(0.1, 0.9, 25)
    result = grouped_ap_contrast(rows, scores, scores, np.ones(25), "double_cold", resamples=30, seed=1)
    assert result["interval"] == [0., 0.]
    assert not result["improvement_supported"]
    small = grouped_ap_contrast(rows.assign(phage_group="one"), scores, scores, np.ones(25),
                                "cold_phage", resamples=30, seed=1)
    assert small["interval"] is None
    assert small["status"] == "descriptive_insufficient_groups"


def test_test_labels_cannot_change_selected_models():
    import runpy
    from pathlib import Path
    fit_scenario = runpy.run_path(str(Path(__file__).resolve().parents[1] /
                                     "scripts/run_paired_baseline.py"))["fit_scenario"]
    hosts, phages, train = fixtures()
    phages.loc["P3"] = [0.5, 1.5]
    phages.loc["P4"] = [1.5, 0.5]
    train = train.assign(partition="fit", regime="fit")
    validation = train.iloc[:2].assign(bacterium=["H1", "H2"], phage="P3",
                                     infection=[0, 1], partition="validation", regime="cold_phage")
    test = validation.assign(phage="P4", partition="test")
    rows = pd.concat([train, validation, test], ignore_index=True).assign(weight=1.)
    rows["host_group"], rows["phage_group"] = rows.bacterium, rows.phage
    plan = {"model_selection": {"logistic_C": [0.01, 0.1], "logistic_max_iter": 2000},
            "evaluation": {"bootstrap_resamples": 5, "bootstrap_seed": 1, "minimum_resampled_groups_for_interval": 5}}
    args = (hosts, phages, ["d"], ["surface"], {p: p for p in phages.index}, plan, "observation_equal")
    original = fit_scenario(rows, *args)[0]
    changed = rows.copy()
    changed.loc[changed.partition == "test", "infection"] = [1, 0]
    assert original == fit_scenario(changed, *args)[0]


def test_nonconvergence_is_a_failure_not_an_alternate_model(monkeypatch):
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    def unconverged(*args, **kwargs):
        warnings.warn("fixture did not converge", ConvergenceWarning)
    monkeypatch.setattr(LogisticRegression, "fit", unconverged)
    with pytest.raises(ConvergenceWarning):
        fit_logistic(np.ones((4, 2)), np.array([0, 1, 0, 1]), np.ones(4), C=.1, max_iter=1)
