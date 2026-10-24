"""Strict dataset registry and provenance verification."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import Any

from .analysis import summarize_interactions
from .data import load_interactions
from .artifacts import sha256_file

MANIFEST_SCHEMA_VERSION = 1
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_GIT_COMMIT = re.compile(r"^[0-9a-f]{40}$")
_TOP_LEVEL_FIELDS = {
    "schema_version",
    "dataset_id",
    "version",
    "source",
    "local_path",
    "sha256",
    "license_status",
    "column_mapping",
    "label_semantics",
    "expected_summary",
}
_SOURCE_FIELDS = {
    "project",
    "repository_url",
    "revision",
    "source_path",
    "download_url",
    "citation",
}
_LABEL_FIELDS = {
    "source_column",
    "positive_value",
    "positive_meaning",
    "negative_value",
    "negative_meaning",
    "operational_definition_status",
    "operational_definition",
    "missing_pair_semantics",
    "model_target_mapping_status",
    "model_target_mapping_note",
}
_SUMMARY_FIELDS = {"interactions", "bacteria", "phages", "positive", "negative"}


@dataclass(frozen=True)
class DatasetVerification:
    """Result of verifying a registered dataset against its manifest."""

    manifest: dict[str, Any]
    path: Path
    sha256: str
    summary: dict[str, Any]


def _require_exact_fields(value: Any, fields: set[str], name: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != fields:
        raise ValueError(f"{name} must contain exactly: {', '.join(sorted(fields))}.")
    return value


def source_endpoint_from_manifest(path: str | Path) -> dict:
    """Describe the source label without asserting an unreviewed assay mapping."""
    manifest = load_dataset_manifest(path)
    labels = manifest["label_semantics"]
    if (labels["positive_value"], labels["negative_value"]) != (1, 0):
        raise ValueError("source endpoint v1 requires explicit binary 1/0 coding")
    return {
        "endpoint_id": f"{manifest['dataset_id']}:source-interaction:v1",
        "kind": "source_defined_binary",
        "definition": labels["operational_definition"],
        "positive_value": 1,
        "negative_value": 0,
        "assay_definition_status": (
            "verified" if labels["operational_definition_status"] == "verified" else "unverified"
        ),
        "source": {
            "dataset_id": manifest["dataset_id"],
            "version": manifest["version"],
            "sha256": manifest["sha256"],
            "label_column": labels["source_column"],
        },
        "evidence": [manifest["source"]["download_url"]],
        "limitations": [
            labels["model_target_mapping_note"],
            "Missing pairs are unknown, not negative; mixture outcomes are outside this endpoint.",
        ],
    }


def _require_text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string.")
    return value


def load_dataset_manifest(path: str | Path) -> dict[str, Any]:
    """Load and strictly validate a versioned dataset manifest."""
    manifest_path = Path(path)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid dataset manifest JSON: {manifest_path}") from exc

    _require_exact_fields(manifest, _TOP_LEVEL_FIELDS, "manifest")
    if manifest["schema_version"] != MANIFEST_SCHEMA_VERSION:
        raise ValueError(
            f"unsupported dataset manifest schema_version: {manifest['schema_version']!r}"
        )
    for field in ("dataset_id", "version", "license_status"):
        _require_text(manifest[field], field)

    source = _require_exact_fields(manifest["source"], _SOURCE_FIELDS, "source")
    for field in _SOURCE_FIELDS:
        _require_text(source[field], f"source.{field}")
    if not _GIT_COMMIT.fullmatch(source["revision"]):
        raise ValueError("source.revision must be a full 40-character Git commit SHA.")

    local_path = Path(_require_text(manifest["local_path"], "local_path"))
    if local_path.is_absolute() or ".." in local_path.parts:
        raise ValueError("local_path must be a repository-relative path without '..'.")

    if not isinstance(manifest["sha256"], str) or not _SHA256.fullmatch(manifest["sha256"]):
        raise ValueError("sha256 must be a lowercase 64-character SHA-256 digest.")

    mapping = manifest["column_mapping"]
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError("column_mapping must be a non-empty object.")
    if any(not isinstance(k, str) or not isinstance(v, str) or not v for k, v in mapping.items()):
        raise ValueError("column_mapping must map non-empty column names to column names.")

    labels = _require_exact_fields(manifest["label_semantics"], _LABEL_FIELDS, "label_semantics")
    for field in _LABEL_FIELDS - {"positive_value", "negative_value"}:
        _require_text(labels[field], f"label_semantics.{field}")
    if type(labels["positive_value"]) is not int or type(labels["negative_value"]) is not int:
        raise ValueError("label values must be integers.")
    if labels["positive_value"] == labels["negative_value"]:
        raise ValueError("positive_value and negative_value must differ.")
    if labels["missing_pair_semantics"] != "unknown_not_negative":
        raise ValueError("missing_pair_semantics must be 'unknown_not_negative'.")
    if labels["model_target_mapping_status"] != "not_verified":
        raise ValueError(
            "model_target_mapping_status must remain 'not_verified' until assay mapping is reviewed."
        )

    summary = _require_exact_fields(manifest["expected_summary"], _SUMMARY_FIELDS, "expected_summary")
    for field, value in summary.items():
        if type(value) is not int or value < 0:
            raise ValueError(f"expected_summary.{field} must be a non-negative integer.")
    if summary["positive"] + summary["negative"] != summary["interactions"]:
        raise ValueError("expected positive and negative counts must sum to interactions.")

    return manifest


def verify_registered_dataset(
    manifest_path: str | Path,
    data_path: str | Path | None = None,
) -> DatasetVerification:
    """Verify file checksum and registered summary statistics."""
    manifest_path = Path(manifest_path)
    manifest = load_dataset_manifest(manifest_path)
    if data_path is None:
        repository_root = manifest_path.resolve().parents[2]
        data_path = repository_root / manifest["local_path"]
    path = Path(data_path)

    if not path.is_file():
        raise FileNotFoundError(f"registered dataset file not found: {path}")

    observed_sha256 = sha256_file(path)
    if observed_sha256 != manifest["sha256"]:
        raise ValueError(
            f"dataset SHA-256 mismatch for {path}: expected {manifest['sha256']}, "
            f"observed {observed_sha256}"
        )

    frame = load_interactions(path)
    summary = summarize_interactions(frame)
    for field, expected in manifest["expected_summary"].items():
        if summary[field] != expected:
            raise ValueError(
                f"dataset summary mismatch for {field}: expected {expected}, "
                f"observed {summary[field]}"
            )

    return DatasetVerification(
        manifest=manifest,
        path=path,
        sha256=observed_sha256,
        summary=summary,
    )
