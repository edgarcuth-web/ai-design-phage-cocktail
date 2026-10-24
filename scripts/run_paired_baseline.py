"""Run the fixed additive paired-model study on every admitted frozen split."""
import argparse
import gzip
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from phage_cocktail.artifacts import atomic_write_json, sha256_file
from phage_cocktail.cohort import canonical_assembly_digest, curate_mlc
from phage_cocktail.data import load_interactions
from phage_cocktail.host_baseline import fit_host_model, predict_host_model
from phage_cocktail.host_study import _features, verify_resources
from phage_cocktail.paired_baseline import (
    FOURMERS, PairEncoder, canonical_fourmers, fit_logistic, predict_logistic,
    sequence_pair_weights, weighted_ap, metric_report, grouped_ap_contrast,
)
from phage_cocktail.registry import verify_registered_dataset

MODELS = {"host_only_pooled_logistic": "host", "phage_only_logistic": "phage",
          "host_phage_additive_logistic": "pair"}


def validate_manifest(rows, observed, host_digests, phage_digests):
    joined = rows.merge(observed, on=["bacterium", "phage"], how="outer", validate="one_to_one", indicator=True)
    if not joined._merge.eq("both").all():
        raise ValueError("frozen manifest must cover exactly the observed cohort")
    joined = joined.drop(columns="_merge")
    fit = joined[joined.partition == "fit"]
    for regime, axes in {"cold_host": ["host_group"], "cold_phage": ["phage_group"],
                         "double_cold": ["host_group", "phage_group"]}.items():
        held = joined[(joined.regime == regime) & joined.partition.isin(["validation", "test"])]
        for axis in axes:
            if set(held[axis]) & set(fit[axis]):
                raise ValueError(f"frozen {regime} has {axis} leakage")
    keys = list(zip(joined.bacterium.map(host_digests), joined.phage.map(phage_digests)))
    if any(None in key for key in keys):
        raise ValueError("missing sequence identity")
    pairs = joined.assign(sequence_pair=pd.Series(keys, index=joined.index))
    if pairs.groupby("sequence_pair")["partition"].nunique().max() > 1:
        raise ValueError("an exact sequence pair crosses partitions")
    return joined


def fit_scenario(rows, hosts, phages, numeric, categorical, phage_digests, plan, weighting):
    fit = rows[rows.partition == "fit"].copy()
    val = rows[(rows.regime == "cold_phage") & (rows.partition == "validation")].copy()
    encoder = PairEncoder.fit(fit, hosts, phages, numeric, categorical, phage_digests)
    selected, trials = {}, []
    for name, kind in MODELS.items():
        x = encoder.transform(fit, hosts, phages, kind)
        vx = encoder.transform(val, hosts, phages, kind)
        best, best_ap = None, -np.inf
        for C in plan["model_selection"]["logistic_C"]:
            model = fit_logistic(x, fit.infection.to_numpy(), fit.weight.to_numpy(),
                                 C=C, max_iter=plan["model_selection"]["logistic_max_iter"])
            ap = weighted_ap(val.infection.to_numpy(), predict_logistic(model, vx), val.weight.to_numpy())
            trials.append(dict(model=name, C=C, validation_cold_phage_ap=ap))
            if ap > best_ap:
                best, best_ap = model, ap
        selected[name] = best
    legacy = fit_host_model(fit, hosts, numeric, categorical, C=0.1, max_iter=2000,
                            sample_weight=None if weighting == "observation_equal" else fit.weight.to_numpy())
    stats = fit.assign(positive_weight=fit.weight * fit.infection).groupby("phage").agg(
        positive=("positive_weight", "sum"), total=("weight", "sum"))
    prevalence = ((stats.positive + 1) / (stats.total + 2)).to_dict()
    global_score = float((fit.weight @ fit.infection + 1) / (fit.weight.sum() + 2))
    artifact = {"schema_version": 1, "kind": "paired_development_models", "encoder": encoder.to_dict(),
                "models": selected, "global_prevalence": global_score, "phage_prevalence": prevalence,
                "legacy_host_logistic": legacy, "weighting": weighting,
                "score_semantics": "uncalibrated_source_mlc_positive"}
    metrics, contrasts, predictions = [], [], []
    for (regime, partition), held in rows[rows.partition.isin(["validation", "test"])].groupby(
            ["regime", "partition"], sort=True):
        scores = {name: predict_logistic(model, encoder.transform(held, hosts, phages, MODELS[name]))
                  for name, model in selected.items()}
        scores["global_prevalence"] = np.full(len(held), global_score)
        scores["phage_prevalence"] = np.array([prevalence.get(p, np.nan) for p in held.phage])
        known = held.phage.isin(legacy["phages"]).to_numpy()
        legacy_scores = np.full(len(held), np.nan)
        if known.any():
            legacy_scores[known] = predict_host_model(legacy, held.loc[known], hosts)
        scores["legacy_host_logistic_refit"] = legacy_scores
        saved = held.copy()
        for name, values in scores.items():
            metrics.append(dict(regime=regime, partition=partition, model=name,
                                **metric_report(held, values, held.weight.to_numpy())))
            saved[name] = values
        predictions.append(saved)
        if partition == "test":
            settings = plan["evaluation"]
            contrast = grouped_ap_contrast(
                held, scores["host_phage_additive_logistic"], scores["host_only_pooled_logistic"],
                held.weight.to_numpy(), regime, resamples=settings["bootstrap_resamples"],
                seed=settings["bootstrap_seed"], minimum_groups=settings["minimum_resampled_groups_for_interval"])
            contrasts.append(dict(regime=regime, **contrast))
    return artifact, trials, metrics, contrasts, pd.concat(predictions, ignore_index=True)


def run(repo, data_root, outdir, *, fetch=False):
    repo, outdir = Path(repo).resolve(), Path(outdir).resolve()
    if outdir.exists():
        raise ValueError("study output directory must be new")
    protocol_path = repo / "data/studies/paired_baseline.v2.json"
    plan = json.loads(protocol_path.read_text())
    parent = plan["parent_v2"]
    if sha256_file(repo / parent["path"]) != parent["sha256"]:
        raise ValueError("parent v2 protocol changed")
    reuse = plan["data_reuse_review"]
    if sha256_file(repo / reuse["path"]) != reuse["sha256"]:
        raise ValueError("data reuse review changed")
    if platform.python_version() != plan["numerical_environment"]["python"]:
        raise ValueError("Python version differs from numerical lock")
    environment = repo / plan["numerical_environment"]["path"]
    if sha256_file(environment) != plan["numerical_environment"]["sha256"]:
        raise ValueError("numerical environment lock changed")
    versions = {"python": platform.python_version()}
    for line in environment.read_text().splitlines():
        if "==" in line and not line.startswith("#"):
            name, version = line.split("==")
            versions[name] = importlib.metadata.version(name)
            if versions[name] != version:
                raise ValueError(f"numerical package version differs: {name}")
    freeze_path = repo / plan["split_freeze"]["report"]
    if sha256_file(freeze_path) != plan["split_freeze"]["sha256"]:
        raise ValueError("split freeze report changed")
    freeze = json.loads(freeze_path.read_text())
    for name, digest in freeze["artifacts_sha256"].items():
        if sha256_file(freeze_path.parent / name) != digest:
            raise ValueError(f"frozen artifact changed: {name}")
    resource_path = repo / plan["inputs"]["resource_lock"]
    if sha256_file(resource_path) != plan["inputs"]["resource_lock_sha256"]:
        raise ValueError("source resource lock changed")
    resources = verify_resources(json.loads(resource_path.read_text()), data_root, fetch=fetch)
    subset = resources[("table", "subset")]
    verify_registered_dataset(repo / plan["inputs"]["subset_manifest"], subset)
    observed, excluded, _ = curate_mlc(load_interactions(subset),
        pd.read_csv(resources[("table", "mlc")], sep=";", index_col=0))
    hosts = sorted(observed.bacterium.unique())
    host_digests = {h: canonical_assembly_digest(resources[("host_genome", h)]) for h in hosts}
    phage_ids = sorted(observed.phage.unique())
    phage_digests = {p: canonical_assembly_digest(resources[("phage_genome", p)]) for p in phage_ids}
    features, numeric = _features(hosts,
        pd.read_csv(resources[("table", "metadata")], sep=";", index_col=0),
        pd.read_csv(resources[("table", "defenses")], sep=";", index_col=0),
        ["O-type", "H-type", "LPS_type", "ABC_serotype"])
    phages = pd.DataFrame([canonical_fourmers(resources[("phage_genome", p)].read_bytes()) for p in phage_ids],
                          index=phage_ids, columns=["4mer:" + s for s in FOURMERS])
    code = [Path(__file__).resolve(), *(repo / "src/phage_cocktail" / name for name in (
        "paired_baseline.py", "host_baseline.py", "cohort.py", "host_study.py", "artifacts.py", "registry.py", "data.py"))]
    provenance = {"protocol_sha256": sha256_file(protocol_path), "versions": versions,
                  "code_base_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
                  "code_sha256": {p.relative_to(repo).as_posix(): sha256_file(p) for p in code},
                  "resource_lock_sha256": sha256_file(resource_path), "split_freeze_report_sha256": sha256_file(freeze_path),
                  "data_reuse_review_sha256": reuse["sha256"]}
    outdir.mkdir(parents=True)
    atomic_write_json(provenance, outdir / "provenance.json")
    phages.to_csv(outdir / "phage_features.csv", index_label="phage")
    features.to_csv(outdir / "host_features.csv", index_label="bacterium")
    excluded.to_csv(outdir / "excluded_pairs.csv", index=False)
    results, metric_rows, trial_rows = [], [], []
    for scenario in freeze["scenarios"]:
        identity, seed = scenario["identity"], scenario["seed"]
        tag = f"i{round(identity*100)}_s{seed}"
        manifest = freeze_path.parent / f"paired_split_{tag}.csv.gz"
        rows = validate_manifest(pd.read_csv(manifest, dtype={"bacterium": str}),
                                 observed, host_digests, phage_digests)
        pair_weights = sequence_pair_weights(rows, host_digests, phage_digests)
        for weighting in ("observation_equal", "sequence_pair_equal"):
            current = rows.assign(weight=np.ones(len(rows)) if weighting == "observation_equal" else pair_weights)
            with threadpool_limits(limits=1):
                artifact, trials, metrics, contrasts, predictions = fit_scenario(
                    current, features, phages, numeric, ["O-type", "H-type", "LPS_type", "ABC_serotype"],
                    phage_digests, plan, weighting)
            scenario_name = f"{tag}_{weighting}"
            artifact.update(scenario=scenario_name, provenance=provenance, split_manifest_sha256=sha256_file(manifest))
            atomic_write_json(artifact, outdir / f"models_{scenario_name}.json")
            payload = predictions.to_csv(index=False, lineterminator="\n").encode()
            (outdir / f"predictions_{scenario_name}.csv.gz").write_bytes(gzip.compress(payload, mtime=0))
            metric_rows.extend(dict(identity=identity, seed=seed, weighting=weighting, **m) for m in metrics)
            trial_rows.extend(dict(identity=identity, seed=seed, weighting=weighting, **t) for t in trials)
            results.append(dict(identity=identity, seed=seed, weighting=weighting,
                                selected_C={name: model["C"] for name, model in artifact["models"].items()},
                                contrasts=contrasts))
            # A checkpoint records completed scenarios; it does not permit fitting a different protocol.
            atomic_write_json({"completed": results}, outdir / "progress.json")
            print(f"Completed {scenario_name}", flush=True)
    pd.DataFrame(metric_rows).to_csv(outdir / "metrics.csv", index=False)
    pd.DataFrame(trial_rows).to_csv(outdir / "selection.csv", index=False)
    summary = {"study_id": plan["study_id"], "provenance": provenance, "scenarios": results,
               "primary_scenario": "i95_s42_observation_equal", "metrics": metric_rows,
               "limitations": ["Source MLC positivity, not verified productive infection or cocktail efficacy.",
                   "Additive compositional pair baseline; no explicit host-phage mechanistic interaction.",
                   "Already inspected development source; no independent release validation, calibration or OOD claim."],
               "artifacts_sha256": {p.name: sha256_file(p) for p in sorted(outdir.iterdir())}}
    atomic_write_json(summary, outdir / "report.json")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    run(Path(__file__).resolve().parents[1], args.data_root, args.outdir, fetch=args.fetch)


if __name__ == "__main__":
    main()
