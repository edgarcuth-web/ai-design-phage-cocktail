"""Pair-benchmark admission must not confuse IDs, genomes or missing labels."""
import pandas as pd
import pytest

from phage_cocktail.pair_cohort import audit_pair_identities


def test_exact_duplicate_phages_keep_labels_and_report_disagreement():
    rows = pd.DataFrame({"bacterium": ["H", "H", "J", "J"],
                         "phage": ["A", "B", "A", "B"], "infection": [0, 1, 1, 1]})
    identities, conflicts = audit_pair_identities(
        rows, {"H": "h", "J": "j"}, {"A": "same", "B": "same"},
    )
    assert identities == {"A": "A", "B": "A"}
    assert len(conflicts) == 2
    assert set(conflicts.bacterium) == {"H"}
    assert set(conflicts.infection) == {0, 1}
    assert len(rows) == 4  # No voting, recoding or row collapse.


def test_identity_coverage_is_exact_and_unknown_pairs_stay_absent():
    rows = pd.DataFrame({"bacterium": ["H", "J"], "phage": ["A", "B"], "infection": [0, 1]})
    identities, conflicts = audit_pair_identities(rows, {"H": "h", "J": "j"}, {"A": "a", "B": "b"})
    assert identities == {"A": "A", "B": "B"}
    assert conflicts.empty
    with pytest.raises(ValueError, match="genome IDs"):
        audit_pair_identities(rows, {"H": "h"}, {"A": "a", "B": "b"})
    with pytest.raises(ValueError, match="genome IDs"):
        audit_pair_identities(rows, {"H": "h", "J": "j"}, {"A": "a", "B": "b", "C": "c"})


def test_identical_host_assemblies_are_one_sequence_pair_identity():
    rows = pd.DataFrame({"bacterium": ["H", "J"], "phage": ["A", "A"], "infection": [0, 1]})
    _, conflicts = audit_pair_identities(rows, {"H": "same", "J": "same"}, {"A": "a"})
    assert len(conflicts) == 2
