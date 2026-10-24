"""Freeze operational phage-neighbour groups and development splits; no fitting."""
import argparse
from hashlib import sha1
import gzip
import json
from pathlib import Path
import subprocess

import pandas as pd

from phage_cocktail.artifacts import (
    atomic_write_json,
    code_provenance,
    read_json_object,
    sha256_file,
    staged_directory,
)
from phage_cocktail.cohort import canonical_assembly_digest, host_aliases, host_key
from phage_cocktail.host_study import verify_resources
from phage_cocktail.pair_cohort import audit_pair_identities
from phage_cocktail.paired_splits import assign_paired_rows, neighbour_groups, parse_paf

# Pinned revision of the original study repository (see paired_source_review.resources.v1.json).
SOURCE_REVISION = "99656d8c020fcd1dbbc30a02ee10274a07655562"
# Primary and sensitivity seeds from data/studies/paired_benchmark.v2.json.
SPLIT_SEEDS = (42, 43, 44)


def review_sources(conflicts, source_files):
    """Trace annotations, never reinterpret them as a new assay endpoint."""
    raw = pd.read_csv(source_files[("table", "raw_interactions")], sep=";", dtype={"bacteria": str})
    metadata = pd.read_csv(source_files[("table", "phage_metadata")], sep=";")
    host_aliases(raw.bacteria.unique())  # Reject ambiguous normalizations.
    if not set(raw.score.unique()) <= {"0", "1", "n"}:
        raise ValueError("unexpected raw score alphabet")
    raw["host_key"] = raw.bacteria.map(host_key)
    selected, comparisons, missing, meta_pairs = [], [], [], []
    for (host, digest), pair in conflicts.groupby(["bacterium", "phage_assembly_sha256"]):
        names = sorted(pair.phage)
        if len(names) != 2:
            raise ValueError("source review expects two-record duplicate phage groups")
        blocks = []
        for phage in names:
            block = raw[(raw.host_key == host_key(host)) & (raw.phage == phage)].copy()
            if block.empty:
                missing.append({"bacterium": host, "phage": phage})
            selected.append(block)
            blocks.append(block)
        contexts = ["host_key", "image", "replicate", "plate", "log_dilution"]
        joined = blocks[0].merge(
            blocks[1], on=contexts, suffixes=("_a", "_b"), validate="one_to_one"
        )
        valid = joined.score_a.isin(["0", "1"]) & joined.score_b.isin(["0", "1"])
        comparisons.append(
            {
                "bacterium": host,
                "phages": names,
                "matched_contexts": len(joined),
                "interpretable_contexts": int(valid.sum()),
                "different_scores": int((valid & joined.score_a.ne(joined.score_b)).sum()),
            }
        )
    for digest, pair in conflicts.groupby("phage_assembly_sha256"):
        names = sorted(pair.phage.unique())
        meta = metadata[metadata.phage.isin(names)]
        if set(meta.phage) != set(names) or meta.phage.duplicated().any():
            raise ValueError("duplicate phages need unique exact metadata matches")
        meta_pairs.append(
            {
                "phages": names,
                "same_recorded_isolation_host": bool(meta.Phage_host.nunique() == 1),
                "same_recorded_species": bool(meta.Species.nunique() == 1),
            }
        )
    traced = pd.concat(selected, ignore_index=True)
    isolation_hosts = set(metadata.Phage_host.map(host_key))
    return {
        "source_revision": SOURCE_REVISION,
        "source_raw_rows": len(raw),
        "conflicting_sequence_pairs": len(comparisons),
        "raw_annotation_rows_traced": len(traced),
        "raw_rows_missing": missing,
        "missing_hosts_listed_as_isolation_strains": sorted(
            {r["bacterium"] for r in missing if host_key(r["bacterium"]) in isolation_hosts}
        ),
        "matched_contexts": sum(r["matched_contexts"] for r in comparisons),
        "interpretable_contexts": sum(r["interpretable_contexts"] for r in comparisons),
        "different_raw_scores_in_matched_context": sum(r["different_scores"] for r in comparisons),
        "pairs_with_different_interpretable_raw_scores": sum(
            r["different_scores"] > 0 for r in comparisons
        ),
        "duplicate_group_metadata": meta_pairs,
        "interpretation": (
            "Different source file IDs/headers and source annotations are preserved. Image "
            "scoring, replicate-to-MLC aggregation and physical stock identity are not "
            "independently verified. No cause or alternate target label is inferred."
        ),
        "missingness": "Raw score n remains uninterpretable, never zero",
        "resources_sha256": {key: sha256_file(path) for (_, key), path in source_files.items()},
    }


def run(repo, preflight, data_root, outdir, minimap2, *, fetch=False):
    repo, preflight, outdir = (Path(path).resolve() for path in (repo, preflight, outdir))
    if outdir.exists():
        raise ValueError("output directory must be new")
    report = read_json_object(preflight / "report.json")
    for name, digest in report["artifacts_sha256"].items():
        path = (preflight / name).resolve()
        if path.parent != preflight or sha256_file(path) != digest:
            raise ValueError("preflight artifact changed or escaped its directory")
    lock_path = repo / "data/studies/genophi_host_baseline.resources.v1.json"
    if sha256_file(lock_path) != report["provenance"]["resource_lock_sha256"]:
        raise ValueError("preflight resource lock differs")
    grouping_path = repo / "data/studies/phage_grouping.v1.json"
    grouping = read_json_object(grouping_path)
    version = subprocess.check_output([str(minimap2), "--version"], text=True).strip()
    if version != grouping["tool"]["version"]:
        raise ValueError("minimap2 version does not match grouping lock")
    lock = read_json_object(lock_path)
    phage_entries = [r for r in lock["resources"] if r["kind"] == "phage_genome"]
    paths = verify_resources({"resources": phage_entries}, data_root, fetch=fetch)
    source_lock_path = repo / "data/studies/paired_source_review.resources.v1.json"
    source_files = verify_resources(read_json_object(source_lock_path), data_root, fetch=fetch)
    id_types = {"bacterium": str, "source_bacterium": str}
    observed = pd.read_csv(preflight / "observed_pairs.csv", dtype=id_types)
    hosts = pd.read_csv(preflight / "host_splits.csv", dtype=id_types)
    host_groups = hosts.set_index("bacterium").host_group.to_dict()
    host_digests = hosts.set_index("bacterium").assembly_sha256.to_dict()
    digests, sequences = {}, {}
    for entry in phage_entries:
        name = entry["id"]
        raw = paths[("phage_genome", name)].read_bytes()
        git_blob = sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if git_blob != entry["git_blob_sha"]:
            raise ValueError("phage bytes differ from native source Git blob")
        lines = raw.decode("ascii").splitlines()
        if sum(line.startswith(">") for line in lines) != 1 or lines[0] != ">" + name:
            raise ValueError("this lock requires single-contig phages with exact source headers")
        sequences[name] = "".join(line.strip() for line in lines[1:]).upper()
        digests[name] = canonical_assembly_digest(paths[("phage_genome", name)])
    _, conflicts = audit_pair_identities(observed, host_digests, digests)
    source_review = review_sources(conflicts, source_files)
    with staged_directory(outdir) as stage:
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
            groups, edges = neighbour_groups(
                digests, blocks, identity=identity, coverage=grouping["minimum_bilateral_coverage"]
            )
            tag = str(round(identity * 100))
            group_table = pd.DataFrame(
                [
                    dict(phage=p, assembly_sha256=digests[p], phage_group=groups[p])
                    for p in sorted(groups)
                ]
            )
            group_table.to_csv(stage / f"phage_groups_i{tag}.csv", index=False)
            edge_table = pd.DataFrame(edges, columns=["left", "right"])
            edge_table.to_csv(stage / f"phage_edges_i{tag}.csv", index=False)
            for seed in SPLIT_SEEDS:
                assigned = assign_paired_rows(
                    observed, host_groups, groups, host_digests, digests, seed=seed
                )
                joined = assigned.merge(
                    observed[["bacterium", "phage", "infection"]],
                    on=["bacterium", "phage"],
                    validate="one_to_one",
                )
                summaries = [
                    dict(
                        regime=regime,
                        partition=part,
                        rows=len(rows),
                        positive=int(rows.infection.sum()),
                        host_groups=int(rows.host_group.nunique()),
                        phage_groups=int(rows.phage_group.nunique()),
                    )
                    for (regime, part), rows in joined.groupby(["regime", "partition"])
                ]
                csv = assigned.to_csv(index=False, lineterminator="\n").encode()
                manifest = stage / f"paired_split_i{tag}_s{seed}.csv.gz"
                manifest.write_bytes(gzip.compress(csv, mtime=0))
                scenarios.append(
                    dict(
                        identity=identity,
                        seed=seed,
                        phage_groups=len(set(groups.values())),
                        edges=len(edges),
                        cells=summaries,
                    )
                )
        atomic_write_json(source_review, stage / "source_review.json")
        fasta_hash = sha256_file(fasta)
        fasta.unlink()  # Do not publish a copy of sequence payloads with the split artifacts.
        compressed_paf = stage / "alignments.paf.gz"
        compressed_paf.write_bytes(gzip.compress(paf.read_bytes(), mtime=0))
        paf.unlink()
        result = {
            "freeze_id": "genophi-paired-splits:v1",
            "training_executed": False,
            "primary_identity": grouping["primary_identity"],
            "primary_seed": SPLIT_SEEDS[0],
            "source_review": source_review,
            "scenarios": scenarios,
            "provenance": {
                "grouping_sha256": sha256_file(grouping_path),
                "preflight_report_sha256": sha256_file(preflight / "report.json"),
                "resource_lock_sha256": sha256_file(lock_path),
                "source_review_lock_sha256": sha256_file(source_lock_path),
                "minimap2_version": version,
                "minimap2_binary_sha256": sha256_file(minimap2),
                "command": ["minimap2", *grouping["arguments"], "phages.fasta", "phages.fasta"],
                "combined_fasta_sha256": fasta_hash,
                **code_provenance(repo, __file__),
            },
            "artifacts_sha256": {
                p.name: sha256_file(p) for p in sorted(stage.iterdir()) if p.name != "minimap2.log"
            },
            "limitations": grouping["limitations"],
        }
        atomic_write_json(result, stage / "report.json")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimap2", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    result = run(
        Path(__file__).resolve().parents[1],
        args.preflight,
        args.data_root,
        args.outdir,
        args.minimap2.resolve(),
        fetch=args.fetch,
    )
    scenarios = [
        {k: s[k] for k in ("identity", "seed", "phage_groups", "edges")} for s in result["scenarios"]
    ]
    print(json.dumps({"source_review": result["source_review"], "scenarios": scenarios}, indent=2))


if __name__ == "__main__":
    main()
