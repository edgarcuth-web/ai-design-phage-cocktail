"""Scientific failure cases for the observed-label, grouped host baseline."""
import json

import numpy as np
import pandas as pd
import pytest

from phage_cocktail.cohort import canonical_assembly_digest, curate_mlc, group_hosts
from phage_cocktail.host_baseline import (
    HostFeatureEncoder, fit_host_model, paired_group_bootstrap, predict_host_model,
)
from phage_cocktail.split import split_by_host_group
from phage_cocktail.host_study import host_features, verify_resources


def test_missing_source_mlc_is_excluded_and_disagreement_is_rejected():
    pairs = pd.DataFrame({"bacterium": ["A1", "B"], "phage": ["P", "P"], "infection": [0, 1]})
    source = pd.DataFrame({"P": [np.nan, 3]}, index=["a-1", "B"])
    observed, excluded, aliases = curate_mlc(pairs, source)
    assert observed.bacterium.tolist() == ["B"]
    assert observed.infection.tolist() == [1]
    assert excluded.bacterium.tolist() == ["A1"]
    assert excluded.reason.tolist() == ["source_mlc_missing"]
    assert aliases == {"A1": "a-1", "B": "B"}
    source.loc["B", "P"] = 0
    with pytest.raises(ValueError, match="disagree"):
        curate_mlc(pairs, source)
    source.index = ["a-1", "A1"]
    with pytest.raises(ValueError, match="ambiguous"):
        curate_mlc(pairs, source)


def test_assembly_identity_ignores_headers_order_and_strand_but_rejects_non_dna(tmp_path):
    left, right = tmp_path / "left.fa", tmp_path / "right.fa"
    left.write_text(">one\nAAACG\n>two\nAGGNN\n")
    right.write_text(">different\nNNCCT\n>also-different\nCGTTT\n")
    assert canonical_assembly_digest(left) == canonical_assembly_digest(right)
    right.write_text(">invalid\nNOT_DNA\n")
    with pytest.raises(ValueError, match="DNA"):
        canonical_assembly_digest(right)


def test_near_host_and_identical_assembly_components_never_cross_partitions():
    hosts = list("ABCDEFGH")
    distances = pd.DataFrame(0.01, index=hosts, columns=hosts)
    np.fill_diagonal(distances.values, 0)
    for a, b in [("A", "B"), ("B", "C")]:
        distances.loc[a, b] = distances.loc[b, a] = 0.00001
    digests = {h: h * 64 for h in hosts}
    digests["D"] = digests["A"]
    groups = group_hosts(hosts, distances, digests, threshold=0.0001)
    assert len({groups[h] for h in "ABCD"}) == 1
    pairs = pd.DataFrame({"bacterium": hosts, "phage": "P", "infection": [0, 1] * 4})
    assigned = split_by_host_group(pairs, groups)
    assert assigned.groupby("host_group").split.nunique().max() == 1
    shuffled = split_by_host_group(pairs.sample(frac=1, random_state=1), groups)
    assert assigned.set_index("bacterium").split.to_dict() == shuffled.set_index("bacterium").split.to_dict()
    with pytest.raises(ValueError, match="digest"):
        group_hosts(hosts, distances, {"A": "a" * 64}, threshold=0.0001)


def test_encoder_uses_training_hosts_only_and_ignores_outcome_metadata():
    train = pd.DataFrame({"d": [0, 2], "O": ["O1", None], "n_infections": [0, 100]}, index=["A", "B"])
    encoder = HostFeatureEncoder.fit(train, ["d"], ["O"])
    before = encoder.to_dict()
    test = pd.DataFrame({"d": [1000], "O": ["TEST_ONLY"], "n_infections": [9999]}, index=["C"])
    result = encoder.transform(test)
    assert before == encoder.to_dict()
    assert "TEST_ONLY" not in json.dumps(before)
    assert result[0, 0] == 999
    assert np.all(result[0, 1:] == 0)
    train.n_infections = -999
    np.testing.assert_array_equal(encoder.transform(train), HostFeatureEncoder.fit(train, ["d"], ["O"]).transform(train))
    with pytest.raises(ValueError, match="finite"):
        encoder.transform(test.assign(d=np.nan))


def test_model_round_trip_and_single_class_fallback_are_explicit():
    features = pd.DataFrame({"d": [0, 1, 2, 3], "O": ["O1", "O1", "O2", "O2"]}, index=list("ABCD"))
    pairs = pd.DataFrame({"bacterium": list("ABCD") * 2, "phage": ["P"] * 4 + ["Q"] * 4, "infection": [0, 0, 1, 1, 0, 0, 0, 0]})
    model = fit_host_model(pairs, features, ["d"], ["O"], C=0.1, max_iter=2000)
    scores = predict_host_model(model, pairs, features)
    np.testing.assert_allclose(scores, predict_host_model(json.loads(json.dumps(model)), pairs, features))
    np.testing.assert_allclose(scores[4:], 1 / 6)
    assert model["phages"]["Q"]["kind"] == "single_class_prevalence"
    with pytest.raises(ValueError, match="unsupported phage"):
        predict_host_model(model, pairs.assign(phage="UNSEEN"), features)


def test_bootstrap_retains_whole_groups_and_has_zero_delta_for_equal_scores():
    pairs = pd.DataFrame({"infection": [0, 1, 0, 1, 1, 0], "host_group": ["A", "A", "B", "B", "C", "C"]})
    scores = np.array([0.1, 0.8, 0.2, 0.7, 0.6, 0.3])
    result = paired_group_bootstrap(pairs, scores, scores, resamples=100, seed=42)
    assert result["delta_average_precision"] == 0
    assert result["interval"] == [0.0, 0.0]
    assert result["groups"] == 3
    assert result["valid_resamples"] == 100


def test_resource_lock_rejects_changed_bytes_and_path_escape(tmp_path):
    from hashlib import sha256
    payload = b"software fixture"
    path = tmp_path / "input.txt"
    path.write_bytes(payload)
    entry = {"id": "input", "kind": "table", "path": "input.txt", "size_bytes": len(payload),
             "sha256": sha256(payload).hexdigest()}
    assert verify_resources({"resources": [entry]}, tmp_path)[("table", "input")] == path
    path.write_bytes(b"different bytes!")
    with pytest.raises(ValueError, match="changed"):
        verify_resources({"resources": [entry]}, tmp_path)
    with pytest.raises(ValueError, match="confined"):
        verify_resources({"resources": [{**entry, "path": "../escape"}]}, tmp_path)
    with pytest.raises(ValueError, match="unique"):
        verify_resources({"resources": [entry, entry]}, tmp_path)


def test_published_feature_join_cannot_include_outcome_columns():
    metadata = pd.DataFrame({"O-type": ["O1"], "n_infections": [93]}, index=["a-1"])
    defenses = pd.DataFrame({"RM": [2]}, index=["A1"])
    features, names = host_features(["A1"], metadata, defenses, ["O-type"])
    assert features.to_dict() == {"O-type": {"A1": "O1"}, "defense:RM": {"A1": 2}}
    assert names == ["defense:RM"]
