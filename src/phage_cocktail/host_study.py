"""Reproducible real-data audit and fixed host-feature development experiment."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import urllib.request

import numpy as np
import pandas as pd
import sklearn

from .artifacts import atomic_write_json, sha256_file
from .benchmark import _metrics, pairwise_baseline_predictions
from .cohort import canonical_assembly_digest, curate_mlc, group_hosts, host_aliases
from .host_baseline import fit_host_model, paired_group_bootstrap, predict_host_model
from .registry import verify_registered_dataset
from .split import split_by_host_group


def verify_resources(lock, data_root, *, fetch=False):
    """Verify a complete file lock; download only when explicitly requested."""
    root = Path(data_root).resolve()
    files, paths, identities = [], set(), set()
    for entry in lock["resources"]:
        relative = Path(entry["path"])
        identity = (entry["kind"], entry["id"])
        if (relative.is_absolute() or ".." in relative.parts or relative.as_posix() in paths
                or identity in identities):
            raise ValueError("resource paths and identities must be confined and unique")
        paths.add(relative.as_posix())
        identities.add(identity)
        path = (root / relative).resolve()
        if not path.is_relative_to(root):
            raise ValueError("resource path escapes the data root")
        files.append((entry, path))

    def admit(item):
        entry, path = item
        if not path.exists() and fetch:
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".partial")
            try:
                digest, size = hashlib.sha256(), 0
                with urllib.request.urlopen(entry["url"], timeout=45) as response, temporary.open("xb") as output:
                    while block := response.read(1024 * 1024):
                        size += len(block)
                        if size > entry["size_bytes"]:
                            raise ValueError(f"resource exceeds locked size: {entry['id']}")
                        digest.update(block)
                        output.write(block)
                if size != entry["size_bytes"] or digest.hexdigest() != entry["sha256"]:
                    raise ValueError(f"download differs from resource lock: {entry['id']}")
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        if (not path.is_file() or path.stat().st_size != entry["size_bytes"]
                or sha256_file(path) != entry["sha256"]):
            raise ValueError(f"missing or changed locked resource: {entry['path']}")
        return ((entry["kind"], entry["id"]), path)

    with ThreadPoolExecutor(max_workers=4) as pool:
        return dict(pool.map(admit, files))


def _features(hosts, metadata, defenses, categorical):
    columns = {}
    for source in (metadata, defenses):
        lookup = host_aliases(source.index)
        source_rows = [lookup[h.replace("-", "").upper()] for h in hosts]
        selected = source.loc[source_rows].copy()
        selected.index = hosts
        if source is metadata:
            for name in categorical:
                columns[name] = selected[name]
        else:
            for name in selected:
                columns[f"defense:{name}"] = selected[name]
    return pd.DataFrame(columns, index=hosts), [f"defense:{name}" for name in defenses]


def _source_provenance(repo, protocol_path, lock_path):
    paths = [Path(__file__), Path(__file__).with_name("cohort.py"),
             Path(__file__).with_name("host_baseline.py"), Path(__file__).with_name("split.py"),
             Path(__file__).with_name("benchmark.py"), Path(__file__).with_name("data.py"),
             Path(__file__).with_name("registry.py"), Path(__file__).with_name("artifacts.py"),
             repo / "scripts/run_host_baseline.py"]
    return {"protocol_sha256": sha256_file(protocol_path), "resource_lock_sha256": sha256_file(lock_path),
            "code_base_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
            "code_files_sha256": {p.relative_to(repo).as_posix(): sha256_file(p) for p in paths},
            "versions": {"python": platform.python_version(), "numpy": np.__version__,
                         "pandas": pd.__version__, "scikit_learn": sklearn.__version__}}


def run_study(repo, data_root, output_dir, *, fetch=False):
    """Freeze the admitted cohort/groups before the first fit or metric calculation."""
    repo, output_dir = Path(repo).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise ValueError("study output directory must be new")
    protocol_path = repo / "data/studies/genophi_host_baseline.v1.json"
    lock_path = repo / "data/studies/genophi_host_baseline.resources.v1.json"
    protocol, lock = json.loads(protocol_path.read_text()), json.loads(lock_path.read_text())
    if protocol["study_id"] != lock["study_id"]:
        raise ValueError("study and resource lock identities differ")
    resources = verify_resources(lock, data_root, fetch=fetch)
    subset = resources[("table", "subset")]
    verify_registered_dataset(repo / protocol["subset_manifest"], subset)
    from .data import load_interactions
    pairs = load_interactions(subset)
    mlc = pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0)
    observed, excluded, aliases = curate_mlc(pairs, mlc)
    hosts = sorted(observed.bacterium.unique())
    if set(hosts) != set(pairs.bacterium):
        raise ValueError("a host has no observed MLC pairs; this protocol requires cohort review")
    assembly_digests = {host: canonical_assembly_digest(resources[("host_genome", host)]) for host in hosts}
    phage_digests = {phage: canonical_assembly_digest(resources[("phage_genome", phage)])
                    for phage in sorted(observed.phage.unique())}
    distances = pd.read_csv(resources[("table", "core_distances")], sep="\t", index_col=0)
    groups = group_hosts(hosts, distances, assembly_digests,
                        threshold=protocol["grouping"]["core_distance_threshold"])
    assigned = split_by_host_group(observed, groups, seed=protocol["split"]["seed"],
                                   fractions=protocol["split"]["fractions"])
    metadata = pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0)
    defenses = pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0)
    categorical = protocol["features"]["categorical"]
    features, numeric = _features(hosts, metadata, defenses, categorical)
    provenance = _source_provenance(repo, protocol_path, lock_path)
    output_dir.mkdir(parents=True)
    assigned.to_csv(output_dir / "observed_pairs.csv", index=False)
    excluded.to_csv(output_dir / "excluded_pairs.csv", index=False)
    host_table = assigned[["bacterium", "host_group", "split"]].drop_duplicates().sort_values("bacterium")
    host_table["source_bacterium"] = host_table.bacterium.map(aliases)
    host_table["assembly_sha256"] = host_table.bacterium.map(assembly_digests)
    host_table.to_csv(output_dir / "host_splits.csv", index=False)
    features.to_csv(output_dir / "host_features.csv", index_label="bacterium")
    atomic_write_json(provenance, output_dir / "provenance.json")
    print(f"Frozen {len(observed)} observed pairs; {len(excluded)} unknown pairs excluded; "
          f"{len(set(groups.values()))} host groups", flush=True)

    train = assigned[assigned.split == "train"].copy()
    for partition in ("validation", "test"):
        held = assigned[assigned.split == partition]
        if set(train.host_group) & set(held.host_group):
            raise ValueError("host group leakage")
        if set(held.phage) - set(train.phage):
            raise ValueError("unsupported held-out phage panel")
    estimator = protocol["estimator"]
    model = fit_host_model(train, features, numeric, categorical,
                           C=estimator["C"], max_iter=estimator["max_iter"])
    model["provenance"] = provenance
    atomic_write_json(model, output_dir / "model.json")
    metrics, prediction_tables = [], []
    for partition in ("validation", "test"):
        held = assigned[assigned.split == partition].copy().reset_index(drop=True)
        predictions = pairwise_baseline_predictions(train, held)
        predictions["host_logistic"] = predict_host_model(model, held, features)
        for name, scores in predictions.items():
            metrics.append({"partition": partition, "model": name, **_metrics(held.infection, scores)})
            held[name] = scores
        prediction_tables.append(held)
    predictions = pd.concat(prediction_tables, ignore_index=True)
    predictions.to_csv(output_dir / "predictions.csv", index=False)
    metrics = pd.DataFrame(metrics)
    metrics.to_csv(output_dir / "metrics.csv", index=False)
    test = predictions[predictions.split == "test"]
    evaluation = protocol["evaluation"]
    interval = paired_group_bootstrap(test, test.host_logistic.to_numpy(), test.phage_prevalence.to_numpy(),
                                      resamples=evaluation["bootstrap_resamples"],
                                      seed=evaluation["bootstrap_seed"],
                                      confidence=evaluation["confidence_level"])
    split_summary = {name: {"hosts": int(part.bacterium.nunique()), "groups": int(part.host_group.nunique()),
                            "pairs": len(part), "positive": int(part.infection.sum())}
                     for name, part in assigned.groupby("split")}
    report = {"study_id": protocol["study_id"], "provenance": provenance,
              "audit": {"candidate_pairs": len(pairs), "observed_pairs": len(observed),
                        "unknown_pairs_excluded": len(excluded), "known_binary_disagreements": 0,
                        "original_missing_encoded_zero": int(excluded.source_binary.eq(0).sum()),
                        "source_mlc_counts": observed.source_mlc.value_counts().sort_index().to_dict(),
                        "missing_by_host": excluded.bacterium.value_counts().to_dict(),
                        "host_genomes_verified": len(assembly_digests), "phage_genomes_verified": len(phage_digests),
                        "unique_host_assemblies": len(set(assembly_digests.values())),
                        "unique_phage_assemblies": len(set(phage_digests.values())),
                        "defense_features": len(numeric), "categorical_features": categorical},
              "split": split_summary, "metrics": metrics.to_dict(orient="records"),
              "primary_comparison": interval,
              "single_class_fallback_phages": [p for p, m in model["phages"].items()
                                                if m["kind"] == "single_class_prevalence"],
              "artifacts_sha256": {p.name: sha256_file(p) for p in sorted(output_dir.iterdir()) if p.is_file()},
              "scope": protocol["output_semantics"]}
    atomic_write_json(report, output_dir / "report.json")
    print(metrics.round(4).to_string(index=False), flush=True)
    print(json.dumps(interval, indent=2), flush=True)
    return report
