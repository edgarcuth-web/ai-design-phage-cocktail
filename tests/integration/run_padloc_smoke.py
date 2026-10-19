"""Run PADLOC v2.0.0 on the upstream PADLOC test genome in CI.

The genome and PADLOC-DB release are pinned by the workflow. This script
constructs a temporary HostState/sequence registry from the downloaded genome,
runs the real external PADLOC executable, and verifies the resulting artifacts.
It does not commit the genome or the annotation output to this repository.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from phage_cocktail.artifacts import canonical_content_sha256
from phage_cocktail.feature_artifact import load_feature_artifact
from phage_cocktail.padloc_adapter import load_padloc_summary, run_padloc_host_defense
from phage_cocktail.sequence_registry import fingerprint_fasta


def write_json(path: Path, document: dict) -> None:
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    fasta = args.fasta.resolve(strict=True)
    input_root = fasta.parent
    asset = fingerprint_fasta(fasta, kind="host_genome")
    now = datetime.now(timezone.utc).isoformat()

    host_state = {
        "schema_version": 1,
        "record_id": "padloc-upstream-smoke-host:v1",
        "host_state_id": "padloc-upstream-smoke-host",
        "revision": 1,
        "host_profile_id": "padloc-upstream-smoke-profile:v1",
        "host_genome_sha256": asset.sequence_sha256,
        "as_of": now,
        "components": [],
        "evidence": [],
        "limitations": [
            "CI integration fixture sourced from the pinned PADLOC repository test data.",
            "This record is for software integration testing, not a curated biological annotation."
        ],
        "previous_state": None,
        "created_at": now,
    }
    host_state["content_sha256"] = canonical_content_sha256(host_state)
    host_state_path = input_root / "host-state.v1.json"
    write_json(host_state_path, host_state)

    sequence_registry = {
        "schema_version": 1,
        "registry_id": "padloc-upstream-smoke-sequence-registry-v1",
        "entries": [
            {
                "kind": "host_genome",
                "sequence_sha256": asset.sequence_sha256,
                "file_sha256": asset.file_sha256,
                "path": fasta.name,
                "status": "available",
            }
        ],
    }
    sequence_registry_path = input_root / "sequence-registry.v1.json"
    write_json(sequence_registry_path, sequence_registry)

    result = run_padloc_host_defense(
        host_state_path=host_state_path,
        sequence_registry_path=sequence_registry_path,
        padloc_data_dir=args.database,
        output_dir=args.outdir,
        cpu=2,
        timeout_seconds=7200,
    )
    summary = load_padloc_summary(result["summary_path"])
    artifact = load_feature_artifact(result["manifest_path"], verify_outputs=True)

    assert summary["subject"]["sequence_sha256"] == asset.sequence_sha256
    assert summary["analyzer"]["version"] == "v2.0.0"
    assert summary["analyzer"]["database"]["version"] == "v2.0.0"
    assert artifact["feature_set"] == "host_defense_systems"
    assert artifact["status"] == "complete"

    print(
        json.dumps(
            {
                "real_padloc_run": True,
                "padloc_version": summary["analyzer"]["version"],
                "padloc_db_version": summary["analyzer"]["database"]["version"],
                "padloc_db_tree_sha256": summary["analyzer"]["database"]["tree_sha256"],
                "outcome": summary["outcome"],
                "systems_detected": summary["systems_detected"],
                "gene_hit_count": summary["gene_hit_count"],
                "feature_artifact_sha256": artifact["content_sha256"],
                "manifest_path": str(result["manifest_path"]),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
