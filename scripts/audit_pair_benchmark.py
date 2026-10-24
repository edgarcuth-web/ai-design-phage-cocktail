"""Re-admit locked paired inputs and emit a no-training benchmark preflight."""
import argparse
from collections import Counter
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd

from phage_cocktail.artifacts import (
    atomic_write_json,
    code_provenance,
    read_json_object,
    sha256_file,
    staged_directory,
)
from phage_cocktail.cohort import canonical_assembly_digest, curate_mlc, group_hosts
from phage_cocktail.data import load_interactions
from phage_cocktail.host_study import HOST_CATEGORICAL_FEATURES, host_features, verify_resources
from phage_cocktail.pair_cohort import audit_pair_identities
from phage_cocktail.registry import verify_registered_dataset
from phage_cocktail.split import split_by_host_group


def run_audit(repo, data_root, outdir, *, fetch=False):
    repo, outdir = Path(repo).resolve(), Path(outdir).resolve()
    if outdir.exists():
        raise ValueError("audit output directory must be new")
    protocol_path = repo / "data/studies/paired_benchmark.v1.json"
    protocol = read_json_object(protocol_path)
    lock_path = repo / protocol["inputs"]["resource_lock"]
    if sha256_file(lock_path) != protocol["inputs"]["resource_lock_sha256"]:
        raise ValueError("paired protocol resource lock changed")
    lock = read_json_object(lock_path)
    resources = verify_resources(lock, data_root, fetch=fetch)
    print(f"Verified {len(resources)} locked resources", flush=True)
    subset = resources[("table", "subset")]
    verify_registered_dataset(repo / protocol["inputs"]["subset_manifest"], subset)
    pairs = load_interactions(subset)
    mlc = pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0)
    observed, excluded, aliases = curate_mlc(pairs, mlc)
    hosts = {
        key: canonical_assembly_digest(path)
        for (kind, key), path in resources.items()
        if kind == "host_genome"
    }
    phages = {
        key: canonical_assembly_digest(path)
        for (kind, key), path in resources.items()
        if kind == "phage_genome"
    }
    phage_groups, conflicts = audit_pair_identities(observed, hosts, phages)
    distances = pd.read_csv(resources[("table", "core_distances")], sep="\t", index_col=0)
    host_groups = group_hosts(
        sorted(hosts),
        distances,
        hosts,
        threshold=protocol["host_grouping"]["core_distance_threshold"],
    )
    assigned = split_by_host_group(observed, host_groups, seed=42)
    host_table = (
        assigned[["bacterium", "host_group", "split"]].drop_duplicates().sort_values("bacterium")
    )
    host_table["source_bacterium"] = host_table.bacterium.map(aliases)
    host_table["assembly_sha256"] = host_table.bacterium.map(hosts)
    archived_path = repo / protocol["inputs"]["archived_host_split"]
    if sha256_file(archived_path) != protocol["inputs"]["archived_host_split_sha256"]:
        raise ValueError("archived host split changed")
    archived = pd.read_csv(archived_path)
    if not host_table.reset_index(drop=True).equals(archived.reset_index(drop=True)):
        raise ValueError("recomputed host identity/grouping differs from archived development split")
    metadata = pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0)
    defenses = pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0)
    features, numeric = host_features(sorted(hosts), metadata, defenses, HOST_CATEGORICAL_FEATURES)
    counts = features[numeric].to_numpy(dtype=float)
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError("invalid defense feature counts")
    summary = {
        "hosts": len(hosts),
        "phages": len(phages),
        "candidate_pairs": len(pairs),
        "observed_pairs": len(observed),
        "unknown_pairs_excluded": len(excluded),
        "unique_host_assemblies": len(set(hosts.values())),
        "unique_phage_assemblies": len(set(phages.values())),
        "host_groups": len(set(host_groups.values())),
        "defense_features": len(numeric),
    }
    for name, expected in protocol["expected_cohort"].items():
        if summary.get(name) != expected:
            raise ValueError(f"cohort changed: {name}: {summary.get(name)} != {expected}")
    print(json.dumps(summary), flush=True)
    phage_table = pd.DataFrame(
        [
            {"phage": name, "assembly_sha256": phages[name], "exact_phage_group": phage_groups[name]}
            for name in sorted(phages)
        ]
    )
    sizes = Counter(phage_groups.values())
    duplicates = {
        group: sorted(p for p, g in phage_groups.items() if g == group)
        for group, size in sizes.items()
        if size > 1
    }
    disagreeing_pairs = conflicts[["host_assembly_sha256", "phage_assembly_sha256"]].drop_duplicates()
    with staged_directory(outdir) as staging:
        assigned.to_csv(staging / "observed_pairs.csv", index=False)
        excluded.to_csv(staging / "excluded_pairs.csv", index=False)
        host_table.to_csv(staging / "host_splits.csv", index=False)
        phage_table.to_csv(staging / "phage_identities.csv", index=False)
        conflicts.to_csv(staging / "sequence_pair_label_disagreements.csv", index=False)
        report = {
            "study_id": protocol["study_id"],
            "protocol_status": protocol["status"],
            "training_executed": False,
            "ready_for_training": False,
            "release_validation_available": False,
            "summary": summary,
            "source_mlc_counts": observed.source_mlc.value_counts().sort_index().to_dict(),
            "duplicate_phage_assemblies": duplicates,
            "sequence_pair_label_disagreement_rows": len(conflicts),
            "sequence_pair_label_disagreement_pairs": len(disagreeing_pairs),
            "disagreement_interpretation": (
                "Preserve all observations; identical genomes do not prove identical stocks "
                "or assay context. No majority vote or recoding."
            ),
            "phage_relatedness_status": "exact_assembly_identity_only",
            "blockers": protocol["before_any_fit"],
            "provenance": {
                "protocol_sha256": sha256_file(protocol_path),
                "resource_lock_sha256": sha256_file(lock_path),
                "subset_manifest_sha256": sha256_file(repo / protocol["inputs"]["subset_manifest"]),
                "archived_host_split_sha256": sha256_file(archived_path),
                **code_provenance(repo, __file__),
                "versions": {
                    "python": platform.python_version(),
                    "numpy": np.__version__,
                    "pandas": pd.__version__,
                },
                "resource_counts_verified": dict(Counter(kind for kind, _ in resources)),
            },
            "artifacts_sha256": {p.name: sha256_file(p) for p in sorted(staging.iterdir())},
        }
        atomic_write_json(report, staging / "report.json")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true", help="Download missing locked public resources")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    report = run_audit(repo, args.data_root, args.outdir, fetch=args.fetch)
    print(f"Preflight complete; training remains blocked: {report['blockers']}")


if __name__ == "__main__":
    main()
