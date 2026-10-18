"""Resolve exact HostState and PhageProfile records from an explicit local registry.

The registry verifies identity and byte/content integrity. It does not establish
biological validity and does not resolve raw genome sequence files or features.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
import json
from pathlib import Path, PurePosixPath

from jsonschema import Draft202012Validator, FormatChecker

from .artifacts import (
    canonical_content_sha256,
    read_json_object,
    read_json_object_with_digest,
)
from .model_interface import validate_model_request

SCHEMA_VERSION = 1
RESOURCE_KINDS = ("host_state", "phage_profile")
REGISTRY_SCHEMA_PATH = (
    Path(__file__).resolve().parent / "schemas" / "resource-registry-v1.schema.json"
)
_HOST_STATE_FIELDS = {
    "content_sha256",
    "schema_version",
    "record_id",
    "host_state_id",
    "revision",
    "host_profile_id",
    "host_genome_sha256",
    "as_of",
    "components",
    "evidence",
    "limitations",
    "previous_state",
    "created_at",
}
_PHAGE_PROFILE_FIELDS = {
    "schema_version",
    "profile_id",
    "phage_id",
    "revision",
    "created_at",
    "provenance",
    "sequence",
    "source_annotations",
}
_SHA256_LENGTH = 64


@dataclass(frozen=True)
class ResourceRegistry:
    """Validated registry manifest and its local resource root."""

    registry_id: str
    path: Path
    root: Path
    file_sha256: str
    entries: tuple[dict, ...]


@dataclass(frozen=True)
class ResolvedModelResources:
    """Exact model input records resolved from a local registry."""

    request_id: str
    request_sha256: str
    registry_id: str
    registry_file_sha256: str
    registry_root: Path
    host_state: dict
    host_state_path: Path
    host_state_file_sha256: str
    phage_profile: dict
    phage_profile_path: Path
    phage_profile_file_sha256: str

    def summary(self) -> dict:
        """Return a portable audit summary without copying biological records."""
        sequence = self.phage_profile.get("sequence")
        phage_sequence_sha256 = (
            sequence.get("sha256")
            if isinstance(sequence, dict)
            else None
        )
        return {
            "schema_version": 1,
            "request_id": self.request_id,
            "request_sha256": self.request_sha256,
            "registry": {
                "registry_id": self.registry_id,
                "file_sha256": self.registry_file_sha256,
            },
            "resources": {
                "host_state": {
                    "resource_id": self.host_state["record_id"],
                    "content_sha256": self.host_state["content_sha256"],
                    "file_sha256": self.host_state_file_sha256,
                    "path": self.host_state_path.relative_to(
                        self.registry_root
                    ).as_posix(),
                },
                "phage_profile": {
                    "resource_id": self.phage_profile["profile_id"],
                    "content_sha256": canonical_content_sha256(
                        self.phage_profile
                    ),
                    "file_sha256": self.phage_profile_file_sha256,
                    "path": self.phage_profile_path.relative_to(
                        self.registry_root
                    ).as_posix(),
                },
            },
            "sequence_references_not_resolved": {
                "host_genome_sha256": self.host_state["host_genome_sha256"],
                "phage_sequence_sha256": phage_sequence_sha256,
            },
            "sequence_assets_resolved": False,
            "biological_validity_established": False,
        }


def _is_sha256(value) -> bool:
    return (
        isinstance(value, str)
        and len(value) == _SHA256_LENGTH
        and all(character in "0123456789abcdef" for character in value)
    )


def _parse_timestamp(value, field_name: str) -> datetime:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.") from None
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{field_name} must include a timezone.")
    return parsed


@lru_cache(maxsize=1)
def _registry_validator() -> Draft202012Validator:
    schema = read_json_object(REGISTRY_SCHEMA_PATH, label="resource registry schema")
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, format_checker=FormatChecker())


def _validate_relative_path(value: str, *, field_name: str) -> PurePosixPath:
    if not isinstance(value, str) or not value.strip() or "\\" in value:
        raise ValueError(f"{field_name} must be a non-empty POSIX relative path.")
    if value.startswith("/"):
        raise ValueError(f"{field_name} must be relative to the registry directory.")
    if any(part in {"", ".", ".."} for part in value.split("/")):
        raise ValueError(
            f"{field_name} must not contain empty, '.' or '..' path components."
        )
    path = PurePosixPath(value)
    if path.is_absolute():
        raise ValueError(f"{field_name} must be relative to the registry directory.")
    return path


def load_resource_registry(path: str | Path) -> ResourceRegistry:
    """Load and strictly validate a local registry manifest.

    Resource paths are always relative to the registry file's directory. Entries
    with the same kind/resource_id or the same path are rejected as ambiguous.
    """
    manifest_path = Path(path)
    try:
        resolved_manifest_path = manifest_path.resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"resource registry not found: {manifest_path}") from exc
    if not resolved_manifest_path.is_file():
        raise ValueError(f"resource registry is not a file: {manifest_path}")

    manifest, registry_file_sha256 = read_json_object_with_digest(
        resolved_manifest_path, label="resource registry"
    )
    if type(manifest.get("schema_version")) is not int or manifest["schema_version"] != 1:
        raise ValueError("resource registry schema_version must be integer 1.")
    if not isinstance(manifest.get("registry_id"), str) or not manifest["registry_id"].strip():
        raise ValueError("resource registry registry_id must be a non-empty string.")
    errors = sorted(
        _registry_validator().iter_errors(manifest),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        summaries = []
        for error in errors[:5]:
            location = ".".join(str(part) for part in error.absolute_path) or "$"
            summaries.append(f"{location}: {error.message}")
        extra = f" (and {len(errors) - 5} more)" if len(errors) > 5 else ""
        raise ValueError(
            "invalid resource registry: " + "; ".join(summaries) + extra
        )

    identifiers = set()
    paths = set()
    entries = []
    for index, entry in enumerate(manifest["entries"]):
        if not entry["resource_id"].strip():
            raise ValueError(f"entries[{index}].resource_id must not be blank.")
        key = (entry["kind"], entry["resource_id"])
        if key in identifiers:
            raise ValueError(
                f"duplicate registry entry for {entry['kind']} "
                f"{entry['resource_id']!r}."
            )
        identifiers.add(key)

        relative_path = _validate_relative_path(
            entry["path"], field_name=f"entries[{index}].path"
        )
        if relative_path in paths:
            raise ValueError(f"duplicate registry path: {entry['path']!r}.")
        paths.add(relative_path)

        if entry["kind"] not in RESOURCE_KINDS:
            raise ValueError(f"unsupported resource kind: {entry['kind']!r}.")
        entries.append(dict(entry))

    return ResourceRegistry(
        registry_id=manifest["registry_id"],
        path=resolved_manifest_path,
        root=resolved_manifest_path.parent,
        file_sha256=registry_file_sha256,
        entries=tuple(entries),
    )


def _resolve_registered_json(
    registry: ResourceRegistry,
    *,
    kind: str,
    resource_id: str,
    expected_content_sha256: str,
) -> tuple[dict, Path, str]:
    matches = [
        entry
        for entry in registry.entries
        if entry["kind"] == kind and entry["resource_id"] == resource_id
    ]
    if not matches:
        raise ValueError(
            f"no registered {kind} resource for exact ID {resource_id!r}."
        )
    entry = matches[0]
    if entry["status"] != "available":
        raise ValueError(
            f"registered {kind} resource {resource_id!r} is "
            f"{entry['status']!r}, not available."
        )
    if entry["content_sha256"] != expected_content_sha256:
        raise ValueError(
            f"registered {kind} resource {resource_id!r} has a different "
            "content SHA-256 from the request."
        )

    relative_path = _validate_relative_path(
        entry["path"], field_name=f"{kind} resource path"
    )
    try:
        resolved_path = (registry.root / Path(*relative_path.parts)).resolve(
            strict=True
        )
    except OSError as exc:
        raise ValueError(
            f"registered {kind} resource file is missing: {entry['path']}"
        ) from exc
    if not resolved_path.is_relative_to(registry.root):
        raise ValueError(
            f"registered {kind} resource path escapes the registry directory."
        )
    if not resolved_path.is_file():
        raise ValueError(
            f"registered {kind} resource path is not a file: {entry['path']}"
        )

    document, observed_file_sha256 = read_json_object_with_digest(
        resolved_path, label=kind
    )
    if observed_file_sha256 != entry["file_sha256"]:
        raise ValueError(
            f"{kind} resource file SHA-256 mismatch for {entry['path']!r}."
        )

    observed_content_sha256 = canonical_content_sha256(document)
    if observed_content_sha256 != entry["content_sha256"]:
        raise ValueError(
            f"{kind} canonical content SHA-256 mismatch for {entry['path']!r}."
        )
    embedded_hash = document.get("content_sha256")
    if embedded_hash is not None and embedded_hash != observed_content_sha256:
        raise ValueError(
            f"{kind} embedded content_sha256 does not match its canonical content."
        )
    return document, resolved_path, observed_file_sha256


def _validate_host_state_record(document: dict, reference: dict) -> None:
    if set(document) != _HOST_STATE_FIELDS:
        raise ValueError("resolved HostState has missing or unknown top-level fields.")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise ValueError("resolved HostState schema_version is not supported.")
    if type(document["revision"]) is not int or document["revision"] < 1:
        raise ValueError("resolved HostState revision must be a positive integer.")
    if document["record_id"] != f"{document['host_state_id']}:v{document['revision']}":
        raise ValueError("resolved HostState record_id does not match its revision.")
    for field in ("components", "evidence", "limitations"):
        if not isinstance(document[field], list):
            raise ValueError(f"resolved HostState {field} must be an array.")
    if document["previous_state"] is not None and not isinstance(
        document["previous_state"], dict
    ):
        raise ValueError("resolved HostState previous_state must be an object or null.")
    _parse_timestamp(document["as_of"], "HostState as_of")
    if _parse_timestamp(document["created_at"], "HostState created_at") < _parse_timestamp(
        document["as_of"], "HostState as_of"
    ):
        raise ValueError("resolved HostState created_at predates as_of.")

    expected_fields = (
        "record_id",
        "content_sha256",
        "host_profile_id",
        "host_genome_sha256",
        "as_of",
    )
    for field in expected_fields:
        if document[field] != reference[field]:
            raise ValueError(
                f"resolved HostState {field} does not match ModelRequest."
            )


def _validate_phage_profile_record(document: dict, request: dict) -> None:
    if set(document) != _PHAGE_PROFILE_FIELDS:
        raise ValueError(
            "resolved PhageProfile has missing or unknown top-level fields."
        )
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise ValueError("resolved PhageProfile schema_version is not supported.")
    if type(document["revision"]) is not int or document["revision"] < 1:
        raise ValueError("resolved PhageProfile revision must be a positive integer.")
    if document["profile_id"] != f"{document['phage_id']}:v{document['revision']}":
        raise ValueError("resolved PhageProfile profile_id does not match its revision.")
    _parse_timestamp(document["created_at"], "PhageProfile created_at")
    if document["profile_id"] != request["phage_profile_id"]:
        raise ValueError("resolved PhageProfile ID does not match ModelRequest.")
    if canonical_content_sha256(document) != request["phage_profile_sha256"]:
        raise ValueError("resolved PhageProfile fingerprint does not match ModelRequest.")

    if not isinstance(document["provenance"], dict):
        raise ValueError("resolved PhageProfile provenance must be an object.")
    if document["provenance"].get("kind") not in {
        "retrieved",
        "natural",
        "engineered",
    }:
        raise ValueError("resolved PhageProfile has unsupported provenance kind.")
    if document["sequence"] is not None:
        sequence = document["sequence"]
        if not isinstance(sequence, dict):
            raise ValueError("resolved PhageProfile sequence must be an object or null.")
        if set(sequence) - {"accession", "sha256", "length"}:
            raise ValueError("resolved PhageProfile sequence has unknown fields.")
        if not sequence.get("accession") and not sequence.get("sha256"):
            raise ValueError(
                "resolved PhageProfile sequence needs an accession or SHA-256."
            )
        if "sha256" in sequence and not _is_sha256(sequence["sha256"]):
            raise ValueError("resolved PhageProfile sequence SHA-256 is invalid.")
        if "length" in sequence and (
            type(sequence["length"]) is not int or sequence["length"] < 1
        ):
            raise ValueError("resolved PhageProfile sequence length is invalid.")
    if not isinstance(document["source_annotations"], dict):
        raise ValueError("resolved PhageProfile source_annotations must be an object.")
    for key, value in document["source_annotations"].items():
        if not isinstance(key, str) or not key.strip():
            raise ValueError("resolved PhageProfile annotation keys must not be blank.")
        if value is not None and not isinstance(value, (str, int, float, bool)):
            raise ValueError("resolved PhageProfile annotations must be scalar values or null.")


def resolve_request_records(
    request: dict,
    registry_path: str | Path,
) -> ResolvedModelResources:
    """Resolve exact HostState and PhageProfile records for a validated request.

    This verifies local registration, exact IDs, canonical content fingerprints,
    raw file digests and request lineage. It does not fetch remote resources,
    resolve genome sequence files, derive features, or assert biological validity.
    """
    request = validate_model_request(request)
    registry = load_resource_registry(registry_path)
    host_reference = request["host_state"]

    host_state, host_path, host_file_sha256 = _resolve_registered_json(
        registry,
        kind="host_state",
        resource_id=host_reference["record_id"],
        expected_content_sha256=host_reference["content_sha256"],
    )
    _validate_host_state_record(host_state, host_reference)

    phage_profile, phage_path, phage_file_sha256 = _resolve_registered_json(
        registry,
        kind="phage_profile",
        resource_id=request["phage_profile_id"],
        expected_content_sha256=request["phage_profile_sha256"],
    )
    _validate_phage_profile_record(phage_profile, request)

    return ResolvedModelResources(
        request_id=request["request_id"],
        request_sha256=request["content_sha256"],
        registry_id=registry.registry_id,
        registry_file_sha256=registry.file_sha256,
        registry_root=registry.root,
        host_state=host_state,
        host_state_path=host_path,
        host_state_file_sha256=host_file_sha256,
        phage_profile=phage_profile,
        phage_profile_path=phage_path,
        phage_profile_file_sha256=phage_file_sha256,
    )
