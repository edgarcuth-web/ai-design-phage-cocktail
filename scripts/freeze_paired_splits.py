"""Freeze operational phage-neighbour groups and development splits; no fitting."""
import argparse
from collections import Counter
from hashlib import sha1
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import pandas as pd

from phage_cocktail.artifacts import atomic_write_json, sha256_file
from phage_cocktail.cohort import canonical_assembly_digest, host_aliases
from phage_cocktail.host_study import verify_resources
from phage_cocktail.pair_cohort import audit_pair_identities
from phage_cocktail.paired_splits import neighbour_groups, assign_paired_rows, parse_paf


def review_sources(conflicts, source_files):
    """Trace annotations, never reinterpret them as a new assay endpoint."""
    raw = pd.read_csv(source_files[("table", "raw_interactions")], sep=";", dtype={"bacteria": str})
    metadata = pd.read_csv(source_files[("table", "phage_metadata")], sep=";")
    host_aliases(raw.bacteria.unique())  # Reject ambiguous normalizations.
    if not set(raw.score.unique()) <= {"0", "1", "n"}:
        raise ValueError("unexpected raw score alphabet")
    normal = lambda value: str(value).replace("-", "").upper()
    raw["host_key"] = raw.bacteria.map(normal)
    selected, comparisons, missing = [], [], []
    meta_pairs = []
    for (host, digest), pair in conflicts.groupby(["bacterium", "phage_assembly_sha256"]):
        names = sorted(pair.phage)
        if len(names) != 2:
            raise ValueError("source review expects the five admitted two-record duplicate groups")
        blocks = []
        for phage in names:
            block = raw[(raw.host_key == normal(host)) & (raw.phage == phage)].copy()
            if block.empty:
                missing.append({"bacterium": host, "phage": phage})
            selected.append(block)
            blocks.append(block)
        contexts = ["host_key", "image", "replicate", "plate", "log_dilution"]
        joined = blocks[0].merge(blocks[1], on=contexts, suffixes=("_a", "_b"), validate="one_to_one")
        valid = joined.score_a.isin(["0", "1"]) & joined.score_b.isin(["0", "1"])
        comparisons.append({"bacterium": host, "phages": names,
                            "matched_contexts": len(joined), "interpretable_contexts": int(valid.sum()),
                            "different_scores": int((valid & joined.score_a.ne(joined.score_b)).sum())})
    for digest, pair in conflicts.groupby("phage_assembly_sha256"):
        names = sorted(pair.phage.unique())
        meta = metadata[metadata.phage.isin(names)]
        if set(meta.phage) != set(names) or meta.phage.duplicated().any():
            raise ValueError("duplicate phages need unique exact metadata matches")
        meta_pairs.append({"phages": names, "same_recorded_isolation_host": bool(meta.Phage_host.nunique() == 1),
                           "same_recorded_species": bool(meta.Species.nunique() == 1)})
    traced = pd.concat(selected, ignore_index=True)
    return {
        "source_revision": "99656d8c020fcd1dbbc30a02ee10274a07655562",
        "source_raw_rows": len(raw), "conflicting_sequence_pairs": len(comparisons),
        "raw_annotation_rows_traced": len(traced), "raw_rows_missing": missing,
        "missing_hosts_listed_as_isolation_strains": sorted({r["bacterium"] for r in missing
            if normal(r["bacterium"]) in set(metadata.Phage_host.map(normal))}),
        "matched_contexts": sum(r["matched_contexts"] for r in comparisons),
        "interpretable_contexts": sum(r["interpretable_contexts"] for r in comparisons),
        "different_raw_scores_in_matched_context": sum(r["different_scores"] for r in comparisons),
        "pairs_with_different_interpretable_raw_scores": sum(r["different_scores"] > 0 for r in comparisons),
        "duplicate_group_metadata": meta_pairs,
        "interpretation": "Different source file IDs/headers and source annotations are preserved. Image scoring, replicate-to-MLC aggregation and physical stock identity are not independently verified. No cause or alternate target label is inferred.",
        "missingness": "Raw score n remains uninterpretable, never zero",
        "resources_sha256": {key: sha256_file(path) for (_, key), path in source_files.items()},
    }


def run(repo, preflight, data_root, outdir, minimap2, *, fetch=False):
    repo, preflight, outdir = Path(repo).resolve(), Path(preflight).resolve(), Path(outdir).resolve()
    if outdir.exists():
        raise ValueError("output directory must be new")
    report = json.loads((preflight / "report.json").read_text())
    for name, digest in report["artifacts_sha256"].items():
        path = (preflight / name).resolve()
        if path.parent != preflight or sha256_file(path) != digest:
            raise ValueError("preflight artifact changed or escaped its directory")
    lock_path = repo / "data/studies/genophi_host_baseline.resources.v1.json"
    if sha256_file(lock_path) != report["provenance"]["resource_lock_sha256"]:
        raise ValueError("preflight resource lock differs")
    grouping_path = repo / "data/studies/phage_grouping.v1.json"
    grouping = json.loads(grouping_path.read_text())
    version = subprocess.check_output([str(minimap2), "--version"], text=True).strip()
    if version != grouping["tool"]["version"]:
        raise ValueError("minimap2 version does not match grouping lock")
    lock = json.loads(lock_path.read_text())
    phage_entries = [r for r in lock["resources"] if r["kind"] == "phage_genome"]
    paths = verify_resources({"resources": phage_entries}, data_root, fetch=fetch)
    source_lock_path = repo / "data/studies/paired_source_review.resources.v1.json"
    source_files = verify_resources(json.loads(source_lock_path.read_text()), data_root, fetch=fetch)
    observed = pd.read_csv(preflight / "observed_pairs.csv", dtype={"bacterium": str, "source_bacterium": str})
    hosts = pd.read_csv(preflight / "host_splits.csv", dtype={"bacterium": str, "source_bacterium": str})
    host_groups = hosts.set_index("bacterium").host_group.to_dict()
    host_digests = hosts.set_index("bacterium").assembly_sha256.to_dict()
    digests, sequences = {}, {}
    for entry in phage_entries:
        name = entry["id"]
        raw = paths[("phage_genome", name)].read_bytes()
        if sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() != entry["git_blob_sha"]:
            raise ValueError("phage bytes differ from native source Git blob")
        lines = raw.decode("ascii").splitlines()
        if sum(line.startswith(">") for line in lines) != 1 or lines[0] != ">" + name:
            raise ValueError("this lock requires single-contig phages with exact source headers")
        sequences[name] = "".join(line.strip() for line in lines[1:]).upper()
        digests[name] = canonical_assembly_digest(paths[("phage_genome", name)])
    _, conflicts = audit_pair_identities(observed, host_digests, digests)
    source_review = review_sources(conflicts, source_files)
    outdir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".freeze-pairs-", dir=outdir.parent))
    try:
        fasta = stage / "phages.fasta"
        fasta.write_text("".join(f">{name}\n{sequences[name]}\n" for name in sorted(sequences)))
        paf, log = stage / "alignments.paf", stage / "minimap2.log"
        command = [str(minimap2), *grouping["arguments"], str(fasta), str(fasta)]
        with paf.open("w") as output, log.open("w") as errors:
            subprocess.run(command, stdout=output, stderr=errors, check=True)
        paf.write_text("\n".join(sorted(paf.read_text().splitlines())) + "\n")
        blocks = parse_paf(paf, {name: len(seq) for name, seq in sequences.items()})
        scenarios = []
        for identity in [grouping["primary_identity"], *grouping["sensitivity_identities"]]:
            groups, edges = neighbour_groups(digests, blocks, identity=identity,
                                             coverage=grouping["minimum_bilateral_coverage"])
            tag = str(round(identity * 100))
            pd.DataFrame([dict(phage=p, assembly_sha256=digests[p], phage_group=groups[p])
                          for p in sorted(groups)]).to_csv(stage / f"phage_groups_i{tag}.csv", index=False)
            pd.DataFrame(edges, columns=["left", "right"]).to_csv(stage / f"phage_edges_i{tag}.csv", index=False)
            for seed in [42, 43, 44]:
                assigned = assign_paired_rows(observed, host_groups, groups, host_digests, digests, seed=seed)
                joined = assigned.merge(observed[["bacterium", "phage", "infection"]],
                                        on=["bacterium", "phage"], validate="one_to_one")
                summaries = []
                for (regime, part), rows in joined.groupby(["regime", "partition"]):
                    summaries.append(dict(regime=regime, partition=part, rows=len(rows),
                                          positive=int(rows.infection.sum()), host_groups=int(rows.host_group.nunique()),
                                          phage_groups=int(rows.phage_group.nunique())))
                csv = assigned.to_csv(index=False, lineterminator="\n").encode()
                (stage / f"paired_split_i{tag}_s{seed}.csv.gz").write_bytes(gzip.compress(csv, mtime=0))
                scenarios.append(dict(identity=identity, seed=seed, phage_groups=len(set(groups.values())),
                                      edges=len(edges), cells=summaries))
        atomic_write_json(source_review, stage / "source_review.json")
        fasta_hash = sha256_file(fasta)
        fasta.unlink()  # Do not publish a copy of sequence payloads with the split artifacts.
        compressed_paf = stage / "alignments.paf.gz"
        compressed_paf.write_bytes(gzip.compress(paf.read_bytes(), mtime=0))
        paf.unlink()
        result = {
            "freeze_id": "genophi-paired-splits:v1", "training_executed": False,
            "primary_identity": grouping["primary_identity"], "primary_seed": 42,
            "source_review": source_review, "scenarios": scenarios,
            "provenance": {"grouping_sha256": sha256_file(grouping_path),
                           "preflight_report_sha256": sha256_file(preflight / "report.json"),
                           "resource_lock_sha256": sha256_file(lock_path),
                           "source_review_lock_sha256": sha256_file(source_lock_path),
                           "minimap2_version": version, "minimap2_binary_sha256": sha256_file(minimap2),
                           "command": ["minimap2", *grouping["arguments"], "phages.fasta", "phages.fasta"],
                           "combined_fasta_sha256": fasta_hash,
                           "code_sha256": {p.relative_to(repo).as_posix(): sha256_file(p) for p in (
                               Path(__file__).resolve(), repo / "src/phage_cocktail/paired_splits.py",
                               repo / "src/phage_cocktail/split.py", repo / "src/phage_cocktail/pair_cohort.py")}},
            "artifacts_sha256": {p.name: sha256_file(p) for p in sorted(stage.iterdir()) if p.name != "minimap2.log"},
            "limitations": grouping["limitations"],
        }
        atomic_write_json(result, stage / "report.json")
        stage.rename(outdir)
        return result
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimap2", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    result = run(Path(__file__).resolve().parents[1], args.preflight, args.data_root,
                 args.outdir, args.minimap2.resolve(), fetch=args.fetch)
    print(json.dumps({"source_review": result["source_review"],
                      "scenarios": [{k: s[k] for k in ("identity", "seed", "phage_groups", "edges")}
                                    for s in result["scenarios"]]}, indent=2))


if __name__ == "__main__":
    main()
