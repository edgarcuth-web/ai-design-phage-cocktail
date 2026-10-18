"""Provenance and integrity contract for computed biological feature outputs.

This module records and verifies artifacts produced by real analysis tools. It
does not run an analyzer, interpret its biology, or claim experimental validity.
Output paths are relative to the FeatureArtifact manifest and are confined to
that directory.
"""

from __future__ import annotations

from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tempfile

from jsonschema import Draft202012Validator, FormatChecker

from .artifacts import (
    canonical_content_sha256,
    canonical_json,
    read_json_object,
    read_json_object_with_digest,
)

SCHEMA_VERSION = 1
SCHEMA_PATH = (
    Path(__file__).resolve().parent / "schemas" / "feature-artifact-v1.schema.json"
)
MAX_MANIFEST_BYTES = 4_000_000
READ_CHUNK_BYTES = 1024 * 1024
OUTPUT_SPEC_FIELDS = {"role", "path", "format", "media_type", "schema_sha256"}


@lru_cache(maxsize=1)
def _feature_artifact_validator() -> Draft202012Validator:
    schema = read_json_object(SCHEMA_PATH, label="FeatureArtifact schema")
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, format_checker=FormatChecker())


def _parameters_sha256(parameters: dict) -> str:
    return hashlib.sha256(canonical_json(parameters).encode("utf-8")).hexdigest()


def _timestamp(value: str, field_name: str) -> datetime:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.") from None
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{field_name} must include a timezone.")
    return parsed


def _relative_path(value: str, *, field_name: str) -> PurePosixPath:
    if not isinstance(value, str) or not value.strip() or "\\" in value:
        raise ValueError(f"{field_name} must be a non-empty POSIX relative path.")
    has_windows_drive = (
        len(value) >= 2 and value[0].isalpha() and value[1] == ":"
    )
    if value.startswith("/") or has_windows_drive or any(
        part in {"", ".", ".."} for part in value.split("/")
    ):
        raise ValueError(
            f"{field_name} must not be absolute or contain empty, '.' or '..' components."
        )
    path = PurePosixPath(value)
    if path.is_absolute():
        raise ValueError(f"{field_name} must be relative.")
    return path


def _manifest_root(manifest_path: str | Path) -> tuple[Path, Path]:
    target = Path(manifest_path).resolve(strict=False)
    root = target.parent.resolve(strict=True)
    if not root.is_dir():
        raise ValueError(f"FeatureArtifact manifest directory does not exist: {root}")
    return target, root


def _resolve_output_path(
    manifest_path: str | Path,
    output_path: str,
    *,
    field_name: str = "output path",
) -> Path:
    relative = _relative_path(output_path, field_name=field_name)
    target, root = _manifest_root(manifest_path)
    try:
        resolved = (root / Path(*relative.parts)).resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"{field_name} does not exist: {output_path!r}") from exc
    if resolved == target:
        raise ValueError("FeatureArtifact output cannot be the manifest itself.")
    if not resolved.is_relative_to(root):
        raise ValueError(f"{field_name} escapes the FeatureArtifact directory.")
    if not resolved.is_file():
        raise ValueError(f"{field_name} is not a regular file: {output_path!r}")
    return resolved


def _digest_and_size(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as stream:
            while chunk := stream.read(READ_CHUNK_BYTES):
                digest.update(chunk)
                size += len(chunk)
    except OSError as exc:
        raise ValueError(f"cannot read FeatureArtifact output {path}: {exc}") from exc
    return digest.hexdigest(), size


def validate_feature_artifact(document: dict) -> dict:
    """Validate schema, canonical manifest hash, parameter hash and run chronology."""
    if not isinstance(document, dict):
        raise ValueError("FeatureArtifact must be a JSON object.")
    errors = sorted(
        _feature_artifact_validator().iter_errors(document),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        summaries = []
        for error in errors[:5]:
            location = ".".join(str(part) for part in error.absolute_path) or "$"
            summaries.append(f"{location}: {error.message}")
        extra = f" (and {len(errors) - 5} more)" if len(errors) > 5 else ""
        raise ValueError("invalid FeatureArtifact: " + "; ".join(summaries) + extra)

    if type(document.get("schema_version")) is not int or document["schema_version"] != SCHEMA_VERSION:
        raise ValueError("FeatureArtifact schema_version must be integer 1.")
    expected_content_sha256 = canonical_content_sha256(document)
    if document.get("content_sha256") != expected_content_sha256:
        raise ValueError(
            "FeatureArtifact content_sha256 mismatch: "
            f"expected {expected_content_sha256}, got {document.get('content_sha256')!r}."
        )
    for field_name in ("artifact_id", "feature_set"):
        if not document[field_name].strip():
            raise ValueError(f"FeatureArtifact {field_name} must not be blank.")
    subject = document["subject"]
    if not subject["id"].strip():
        raise ValueError("FeatureArtifact subject.id must not be blank.")
    if subject.get("accession") is not None and not subject["accession"].strip():
        raise ValueError("FeatureArtifact subject.accession must not be blank.")
    tool = document["analysis"]["tool"]
    for field_name in ("name", "version"):
        if not tool[field_name].strip():
            raise ValueError(f"FeatureArtifact analysis.tool.{field_name} must not be blank.")
    for reference in document["analysis"]["reference_data"]:
        for field_name in ("name", "version"):
            if not reference[field_name].strip():
                raise ValueError(
                    f"FeatureArtifact reference_data.{field_name} must not be blank."
                )
    for output in document["outputs"]:
        for field_name in ("role", "path", "format", "media_type"):
            if not output[field_name].strip():
                raise ValueError(f"FeatureArtifact output {field_name} must not be blank.")
    if any(not item.strip() for item in document["limitations"]):
        raise ValueError("FeatureArtifact limitations must not contain blank entries.")

    registry_fields = (
        subject.get("sequence_registry_id"),
        subject.get("sequence_registry_file_sha256"),
    )
    if any(value is not None for value in registry_fields) and not all(
        isinstance(value, str) and value.strip() for value in registry_fields
    ):
        raise ValueError(
            "subject sequence_registry_id and sequence_registry_file_sha256 "
            "must be provided together."
        )

    analysis = document["analysis"]
    expected_parameters_sha256 = _parameters_sha256(analysis["parameters"])
    if analysis["parameters_sha256"] != expected_parameters_sha256:
        raise ValueError(
            "FeatureArtifact parameters_sha256 mismatch: "
            f"expected {expected_parameters_sha256}, got "
            f"{analysis['parameters_sha256']!r}."
        )

    started_at = _timestamp(analysis["started_at"], "analysis.started_at")
    finished_at = _timestamp(analysis["finished_at"], "analysis.finished_at")
    created_at = _timestamp(document["created_at"], "created_at")
    if finished_at < started_at:
        raise ValueError("analysis.finished_at cannot predate analysis.started_at.")
    if created_at < finished_at:
        raise ValueError("FeatureArtifact created_at cannot predate analysis.finished_at.")

    output_paths = set()
    for output in document["outputs"]:
        relative = _relative_path(output["path"], field_name="output path")
        if relative in output_paths:
            raise ValueError(f"duplicate FeatureArtifact output path: {output['path']!r}.")
        output_paths.add(relative)

    reference_identities = set()
    for reference in analysis["reference_data"]:
        identity = (reference["name"], reference["version"])
        if identity in reference_identities:
            raise ValueError(f"duplicate reference snapshot entry: {identity!r}.")
        reference_identities.add(identity)

    return document


def verify_feature_artifact_outputs(
    document: dict, manifest_path: str | Path
) -> dict:
    """Verify each output's exact byte digest and size under the manifest directory."""
    validate_feature_artifact(document)
    verified = []
    for output in document["outputs"]:
        path = _resolve_output_path(
            manifest_path, output["path"], field_name="output path"
        )
        observed_sha256, observed_size = _digest_and_size(path)
        if observed_sha256 != output["sha256"]:
            raise ValueError(
                f"FeatureArtifact output SHA-256 mismatch for {output['path']!r}: "
                f"expected {output['sha256']}, got {observed_sha256}."
            )
        if observed_size != output["size_bytes"]:
            raise ValueError(
                f"FeatureArtifact output size mismatch for {output['path']!r}: "
                f"expected {output['size_bytes']}, got {observed_size}."
            )
        verified.append(
            {
                "role": output["role"],
                "path": output["path"],
                "sha256": observed_sha256,
                "size_bytes": observed_size,
            }
        )
    return {"verified": True, "outputs": verified}


def build_feature_artifact(
    *,
    manifest_path: str | Path,
    artifact_id: str,
    feature_set: str,
    subject: dict,
    tool: dict,
    parameters: dict,
    reference_data: list[dict] | tuple[dict, ...],
    outputs: list[dict] | tuple[dict, ...],
    started_at: str,
    finished_at: str,
    limitations: list[str] | tuple[str, ...],
    status: str = "complete",
    environment: dict | None = None,
    created_at: str | None = None,
) -> dict:
    """Build a sealed manifest from existing tool output files.

    Output specs contain role/path/format/media_type and optional schema_sha256.
    This function hashes existing output files; it never creates feature values
    or runs a biological analyzer.
    """
    if not isinstance(parameters, dict):
        raise ValueError("analysis parameters must be a JSON object.")
    analysis = {
        "tool": dict(tool),
        "started_at": started_at,
        "finished_at": finished_at,
        "parameters": parameters,
        "parameters_sha256": _parameters_sha256(parameters),
        "reference_data": [dict(item) for item in reference_data],
    }
    if environment is not None:
        if not isinstance(environment, dict):
            raise ValueError("analysis environment must be a JSON object.")
        analysis["environment"] = dict(environment)

    prepared_outputs = []
    for index, spec in enumerate(outputs):
        if not isinstance(spec, dict):
            raise ValueError(f"outputs[{index}] must be an object.")
        unknown = set(spec) - OUTPUT_SPEC_FIELDS
        missing = {"role", "path", "format", "media_type"} - set(spec)
        if unknown or missing:
            raise ValueError(
                f"outputs[{index}] fields invalid; missing={sorted(missing)}, "
                f"unknown={sorted(unknown)}."
            )
        output = dict(spec)
        path = _resolve_output_path(
            manifest_path, output["path"], field_name=f"outputs[{index}].path"
        )
        output["sha256"], output["size_bytes"] = _digest_and_size(path)
        prepared_outputs.append(output)

    document = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": artifact_id,
        "created_at": created_at or datetime.now(timezone.utc).isoformat(),
        "feature_set": feature_set,
        "subject": dict(subject),
        "analysis": analysis,
        "outputs": prepared_outputs,
        "status": status,
        "limitations": list(limitations),
    }
    document["content_sha256"] = canonical_content_sha256(document)
    return validate_feature_artifact(document)


def load_feature_artifact(
    path: str | Path, *, verify_outputs: bool = True
) -> dict:
    """Load and validate a manifest; optionally re-hash all referenced output files."""
    manifest_path = Path(path)
    try:
        resolved_manifest = manifest_path.resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"FeatureArtifact manifest not found: {manifest_path}") from exc
    if not resolved_manifest.is_file():
        raise ValueError(f"FeatureArtifact manifest is not a file: {manifest_path}")
    if resolved_manifest.stat().st_size > MAX_MANIFEST_BYTES:
        raise ValueError(
            f"FeatureArtifact manifest exceeds the {MAX_MANIFEST_BYTES}-byte limit."
        )
    document, _ = read_json_object_with_digest(
        resolved_manifest, label="FeatureArtifact manifest"
    )
    validate_feature_artifact(document)
    if verify_outputs:
        verify_feature_artifact_outputs(document, resolved_manifest)
    return document


def write_feature_artifact(
    document: dict, path: str | Path, *, verify_outputs: bool = True
) -> None:
    """Atomically write a validated FeatureArtifact manifest."""
    validate_feature_artifact(document)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if verify_outputs:
        verify_feature_artifact_outputs(document, target)
    payload = json.dumps(document, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    temporary_name = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=target.parent,
            prefix=f".{target.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, target)
    except OSError:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except OSError:
                pass
        raise
