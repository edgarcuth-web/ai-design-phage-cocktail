"""Leakage failures that must block a paired-model development benchmark."""
import pandas as pd
import pytest

from phage_cocktail.paired_splits import neighbour_groups, assign_paired_rows, parse_paf


def hit(q, t, start=0, end=900, matches=890, length=900):
    return dict(query=q, target=t, qlen=1000, tlen=1000, qs=start, qe=end,
                ts=start, te=end, matches=matches, block=length, strand="+")


def test_coverage_is_bilateral_and_repeats_cannot_inflate_it():
    digests = {"A": "a", "B": "b", "C": "c", "D": "c"}
    blocks = [hit("A", "B", 0, 500, 500, 500), hit("A", "B", 0, 500, 500, 500)]
    groups, edges = neighbour_groups(digests, blocks, identity=0.95, coverage=0.85)
    assert groups["A"] != groups["B"]
    assert groups["C"] == groups["D"]
    assert edges == [("C", "D")]
    blocks = [{**hit("A", "B"), "tlen": 2000}]
    groups, _ = neighbour_groups(digests, blocks, identity=0.95, coverage=0.85)
    assert groups["A"] != groups["B"]


def test_edges_are_direction_invariant_and_transitive():
    digests = {"A": "a", "B": "b", "C": "c"}
    groups, edges = neighbour_groups(digests, [hit("B", "A"), hit("B", "C")], identity=0.95, coverage=0.85)
    assert groups == {"A": "A", "B": "A", "C": "A"}
    assert edges == [("A", "B"), ("B", "C")]
    assert neighbour_groups(digests, [hit("B", "C"), hit("B", "A")], identity=0.95, coverage=0.85)[0] == groups


def test_unknown_alignment_identity_is_rejected():
    with pytest.raises(ValueError, match="unknown"):
        neighbour_groups({"A": "a"}, [hit("A", "UNKNOWN")], identity=0.95, coverage=0.85)


def test_approximate_paf_or_wrong_genome_lengths_cannot_define_groups(tmp_path):
    paf = tmp_path / "pairs.paf"
    line = "A\t1000\t0\t900\t+\tB\t1000\t0\t900\t890\t900\t60"
    paf.write_text(line + "\n")
    with pytest.raises(ValueError, match="base-level"):
        parse_paf(paf, {"A": 1000, "B": 1000})
    paf.write_text(line + "\tcg:Z:900M\n")
    assert len(parse_paf(paf, {"A": 1000, "B": 1000})) == 1
    with pytest.raises(ValueError, match="lengths"):
        parse_paf(paf, {"A": 1001, "B": 1000})


def cohort():
    # Large enough panel for deterministic masks to retain every intended seen group.
    rows = pd.DataFrame([(f"H{i}", f"P{j}", (i+j)%2) for i in range(20) for j in range(20)],
                        columns=["bacterium", "phage", "infection"])
    hg = {f"H{i}": f"H{i}" for i in range(20)}
    pg = {f"P{i}": f"P{i}" for i in range(20)}
    return rows, hg, pg


def test_all_rectangles_are_disjoint_and_cold_axes_stay_unseen():
    rows, hg, pg = cohort()
    assigned = assign_paired_rows(rows, hg, pg, hg, pg, seed=42)
    fit = assigned[assigned.partition == "fit"]
    assert len(assigned) == len(rows)
    assert not assigned[["bacterium", "phage"]].duplicated().any()
    for regime in ("cold_host", "cold_phage", "double_cold"):
        for part in ("validation", "test"):
            held = assigned[(assigned.regime == regime) & (assigned.partition == part)]
            assert not held.empty
            if regime in ("cold_host", "double_cold"):
                assert not set(held.host_group) & set(fit.host_group)
            if regime in ("cold_phage", "double_cold"):
                assert not set(held.phage_group) & set(fit.phage_group)
    assert set(assigned.partition) == {"fit", "validation", "test", "unused"}
    shuffled = assign_paired_rows(rows.sample(frac=1, random_state=4), hg, pg, hg, pg, seed=42)
    pd.testing.assert_frame_equal(assigned, shuffled)


def test_duplicate_genomes_share_groups_and_warm_masks():
    rows, hg, pg = cohort()
    extra = rows[rows.phage == "P0"].assign(phage="COPY")
    rows = pd.concat([rows, extra], ignore_index=True)
    pg["COPY"] = pg["P0"]
    pdig = {p: p for p in pg}
    pdig["COPY"] = pdig["P0"]
    assigned = assign_paired_rows(rows, hg, pg, hg, pdig, seed=42)
    paired = assigned[assigned.phage.isin(["P0", "COPY"])]
    assert paired.groupby("bacterium").partition.nunique().max() == 1
    assert paired.groupby("bacterium").regime.nunique().max() == 1
    bad_groups = {**pg, "COPY": "different"}
    with pytest.raises(ValueError, match="identical.*phage"):
        assign_paired_rows(rows, hg, bad_groups, hg, pdig, seed=42)


def test_one_class_fit_is_not_silently_repaired():
    rows, hg, pg = cohort()
    with pytest.raises(ValueError, match="both observed classes"):
        assign_paired_rows(rows.assign(infection=0), hg, pg, hg, pg, seed=42)


def test_uninterpretable_raw_annotations_never_become_negative_labels(tmp_path):
    import runpy
    from pathlib import Path
    review = runpy.run_path(str(Path(__file__).resolve().parents[1] /
                               "scripts/freeze_paired_splits.py"))["review_sources"]
    raw_path, meta_path = tmp_path / "raw.csv", tmp_path / "meta.csv"
    pd.DataFrame({"bacteria": ["H", "H", "H", "H"], "phage": ["A", "B", "A", "B"],
                  "image": ["one", "one", "two", "two"], "replicate": [1, 1, 2, 2],
                  "plate": "A", "log_dilution": 0, "score": ["n", "0", "1", "0"]}).to_csv(
        raw_path, sep=";", index=False)
    pd.DataFrame({"phage": ["A", "B"], "Phage_host": ["H", "H"], "Species": ["S", "S"]}).to_csv(
        meta_path, sep=";", index=False)
    conflicts = pd.DataFrame({"bacterium": ["H", "H"], "phage": ["A", "B"],
                              "phage_assembly_sha256": ["same", "same"]})
    result = review(conflicts, {("table", "raw_interactions"): raw_path,
                                ("table", "phage_metadata"): meta_path})
    assert result["matched_contexts"] == 2
    assert result["interpretable_contexts"] == 1
    assert result["different_raw_scores_in_matched_context"] == 1
