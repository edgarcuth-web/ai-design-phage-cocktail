"""Resolve exact host-genome and phage-sequence FASTA assets.

Two digests are deliberately kept separate:
- file_sha256 hashes the exact original FASTA bytes;
- sequence_sha256 hashes normalized contigs using the lab's GenomeSummary
  convention: uppercase IUPAC DNA bases, headers/whitespace omitted, and a NUL
  byte before each contig.

Successful resolution proves local identity and integrity only, not biological
validity or correctness of the source annotations.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
import io
from pathlib import Path, PurePosixPath

from jsonschema import Draft202012Validator, FormatChecker

from .artifacts import read_json_object, read_json_object_with_digest
from .resource_registry import ResolvedModelResources, resolve_request_records

SCHEMA_VERSION = 1
SEQUENCE_KINDS = ("host_genome", "phage_sequence")
MAX_BASES = 20_000_000
MAX_FILE_BYTES = MAX_BASES * 2
DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
SHA256_LENGTH = 64
SCHEMA_PATH = (
    Path(__file__).resolve().parent / "schemas" / "sequence-registry-v1.schema.json"
)


@dataclass(frozen=True)
class SequenceRegistry:
    """Validated local index of content-addressed sequence files."""

    registry_id: str
    path: Path
    root: Path
    file_sha256: str
    entries: tuple[dict, ...]


@dataclass(frozen=True)
class ResolvedSequenceAsset:
    """One verified FASTA asset; sequence contents are not retained."""

    kind: str
    path: Path
    file_sha256: str
    sequence_sha256: str
    bases: int
    contigs: int
    accession: str | None = None

    def summary(self, root: Path) -> dict:
        result = {
            "kind": self.kind,
            "sequence_sha256": self.sequence_sha256,
            "file_sha256": self.file_sha256,
            "bases": self.bases,
            "contigs": self.contigs,
            "path": self.path.relative_to(root).as_posix(),
        }
        if self.accession is not None:
            result["accession"] = self.accession
        return result


@dataclass(frozen=True)
class ResolvedModelSequences:
    """Verified sequence assets tied to exact model input records."""

    resources: ResolvedModelResources
    sequence_registry_id: str
    sequence_registry_file_sha256: str
    sequence_registry_root: Path
    host_genome: ResolvedSequenceAsset
    phage_sequence: ResolvedSequenceAsset

    def summary(self) -> dict:
        return {
            "schema_version": 1,
            "request_id": self.resources.request_id,
            "request_sha256": self.resources.request_sha256,
            "record_resolution": self.resources.summary(),
            "sequence_registry": {
                "registry_id": self.sequence_registry_id,
                "file_sha256": self.sequence_registry_file_sha256,
            },
            "sequences": {
                "host_genome": self.host_genome.summary(self.sequence_registry_root),
                "phage_sequence": self.phage_sequence.summary(
                    self.sequence_registry_root
                ),
            },
            "sequence_assets_resolved": True,
            "features_derived": False,
            "biological_validity_established": False,
        }


def _is_sha256(value) -> bool:
    return (
        isinstance(value, str)
        and len(value) == SHA256_LENGTH
        and all(character in "0123456789abcdef" for character in value)
    )


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


@lru_cache(maxsize=1)
def _sequence_registry_validator() -> Draft202012Validator:
    schema = read_json_object(SCHEMA_PATH, label="sequence registry schema")
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, format_checker=FormatChecker())


def load_sequence_registry(path: str | Path) -> SequenceRegistry:
    """Load and validate a local sequence registry manifest."""
    manifest_path = Path(path)
    try:
        resolved_manifest_path = manifest_path.resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"sequence registry not found: {manifest_path}") from exc
    if not resolved_manifest_path.is_file():
        raise ValueError(f"sequence registry is not a file: {manifest_path}")

    manifest, manifest_sha256 = read_json_object_with_digest(
        resolved_manifest_path, label="sequence registry"
    )
    if (
        type(manifest.get("schema_version")) is not int
        or manifest["schema_version"] != SCHEMA_VERSION
    ):
        raise ValueError("sequence registry schema_version must be integer 1.")
    if (
        not isinstance(manifest.get("registry_id"), str)
        or not manifest["registry_id"].strip()
    ):
        raise ValueError("sequence registry registry_id must be a non-empty string.")

    errors = sorted(
        _sequence_registry_validator().iter_errors(manifest),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        summaries = []
        for error in errors[:5]:
            location = ".".join(str(part) for part in error.absolute_path) or "$"
            summaries.append(f"{location}: {error.message}")
        extra = f" (and {len(errors) - 5} more)" if len(errors) > 5 else ""
        raise ValueError(
            "invalid sequence registry: " + "; ".join(summaries) + extra
        )

    identities = set()
    paths = set()
    accessions = set()
    entries = []
    for index, entry in enumerate(manifest["entries"]):
        key = (entry["kind"], entry["sequence_sha256"])
        if key in identities:
            raise ValueError(
                f"duplicate sequence registry entry for {entry['kind']} "
                f"{entry['sequence_sha256']!r}."
            )
        identities.add(key)

        relative_path = _validate_relative_path(
            entry["path"], field_name=f"entries[{index}].path"
        )
        if relative_path in paths:
            raise ValueError(f"duplicate sequence registry path: {entry['path']!r}.")
        paths.add(relative_path)

        if entry["kind"] == "host_genome" and "accession" in entry:
            raise ValueError("host_genome entries must not declare a phage accession.")
        if entry["kind"] == "phage_sequence" and "accession" in entry:
            accession = entry["accession"]
            if accession in accessions:
                raise ValueError(
                    f"duplicate phage sequence accession {accession!r}."
                )
            accessions.add(accession)
        entries.append(dict(entry))

    return SequenceRegistry(
        registry_id=manifest["registry_id"],
        path=resolved_manifest_path,
        root=resolved_manifest_path.parent,
        file_sha256=manifest_sha256,
        entries=tuple(entries),
    )


def _parse_fasta_bytes(raw: bytes, *, label: str) -> tuple[str, int, int]:
    """Apply the same normalized-contig digest semantics as phage_lab/biology.py."""
    if len(raw) > MAX_FILE_BYTES:
        raise ValueError(f"{label} FASTA exceeds the {MAX_FILE_BYTES}-byte file limit.")
    try:
        decoded = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{label} FASTA must be UTF-8 text: {exc}") from exc

    digest = hashlib.sha256()
    bases = 0
    contigs = 0
    current_bases = 0

    # newline=None mirrors Python text-mode universal-newline handling used by
    # the lab's read_fasta implementation.
    for line_number, raw_line in enumerate(io.StringIO(decoded, newline=None), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError(f"{label}: empty FASTA header at line {line_number}.")
            if contigs and current_bases == 0:
                raise ValueError(f"{label}: FASTA contains an empty contig.")
            digest.update(b"\0")
            contigs += 1
            current_bases = 0
            continue
        if not contigs:
            raise ValueError(f"{label}: FASTA sequence must follow a >header line.")

        sequence = line.upper()
        if set(sequence) - DNA_ALPHABET:
            raise ValueError(
                f"{label}: invalid DNA characters at line {line_number}."
            )
        bases += len(sequence)
        current_bases += len(sequence)
        if bases > MAX_BASES:
            raise ValueError(f"{label}: FASTA exceeds the {MAX_BASES}-base limit.")
        digest.update(sequence.encode("ascii"))

    if not contigs or not current_bases:
        raise ValueError(f"{label}: FASTA has no sequence or ends with an empty contig.")
    return digest.hexdigest(), bases, contigs


def _read_sequence_file(path: Path, *, label: str) -> tuple[str, str, int, int]:
    """Hash and parse the same bounded byte buffer to avoid read/hash races."""
    try:
        with path.open("rb") as stream:
            raw = stream.read(MAX_FILE_BYTES + 1)
    except OSError as exc:
        raise ValueError(f"cannot read {label} sequence file {path}: {exc}") from exc
    if len(raw) > MAX_FILE_BYTES:
        raise ValueError(f"{label} FASTA exceeds the {MAX_FILE_BYTES}-byte file limit.")
    file_sha256 = hashlib.sha256(raw).hexdigest()
    sequence_sha256, bases, contigs = _parse_fasta_bytes(raw, label=label)
    return file_sha256, sequence_sha256, bases, contigs


def _resolve_sequence(
    registry: SequenceRegistry,
    *,
    kind: str,
    expected_sequence_sha256: str,
    expected_accession: str | None = None,
    expected_length: int | None = None,
) -> ResolvedSequenceAsset:
    if not _is_sha256(expected_sequence_sha256):
        raise ValueError(f"{kind} requires an explicit normalized sequence SHA-256.")

    matches = [
        entry
        for entry in registry.entries
        if entry["kind"] == kind
        and entry["sequence_sha256"] == expected_sequence_sha256
    ]
    if not matches:
        raise ValueError(
            f"no registered {kind} sequence for exact normalized SHA-256 "
            f"{expected_sequence_sha256!r}."
        )
    entry = matches[0]
    if entry["status"] != "available":
        raise ValueError(
            f"registered {kind} sequence {expected_sequence_sha256!r} is "
            f"{entry['status']!r}, not available."
        )
    if expected_accession is not None and entry.get("accession") != expected_accession:
        raise ValueError(
            f"registered {kind} sequence accession does not match the profile reference."
        )

    relative_path = _validate_relative_path(
        entry["path"], field_name=f"{kind} sequence path"
    )
    try:
        resolved_path = (registry.root / Path(*relative_path.parts)).resolve(
            strict=True
        )
    except OSError as exc:
        raise ValueError(
            f"registered {kind} sequence file is missing: {entry['path']}"
        ) from exc
    if not resolved_path.is_relative_to(registry.root):
        raise ValueError(
            f"{kind} sequence path escapes the sequence registry directory."
        )
    if not resolved_path.is_file():
        raise ValueError(
            f"registered {kind} sequence path is not a file: {entry['path']}"
        )

    file_sha256, observed_sequence_sha256, bases, contigs = _read_sequence_file(
        resolved_path, label=kind
    )
    if file_sha256 != entry["file_sha256"]:
        raise ValueError(
            f"{kind} sequence file SHA-256 mismatch for {entry['path']!r}."
        )
    if observed_sequence_sha256 != entry["sequence_sha256"]:
        raise ValueError(
            f"{kind} normalized sequence SHA-256 mismatch for {entry['path']!r}."
        )
    if observed_sequence_sha256 != expected_sequence_sha256:
        raise ValueError(f"{kind} sequence does not match the requested fingerprint.")
    if expected_length is not None and bases != expected_length:
        raise ValueError(
            f"{kind} sequence length mismatch: profile says {expected_length}, "
            f"observed {bases} bases."
        )

    return ResolvedSequenceAsset(
        kind=kind,
        path=resolved_path,
        file_sha256=file_sha256,
        sequence_sha256=observed_sequence_sha256,
        bases=bases,
        contigs=contigs,
        accession=entry.get("accession"),
    )


def resolve_request_sequences(
    request: dict,
    resource_registry_path: str | Path,
    sequence_registry_path: str | Path,
) -> ResolvedModelSequences:
    """Resolve both sequence assets after validating exact input records.

    Accession-only PhageProfile sequence references fail closed: without an
    expected sequence SHA-256, the intended sequence bytes cannot be verified.
    No remote fetch, feature extraction, inference, or biological validation is
    performed.
    """
    resources = resolve_request_records(request, resource_registry_path)
    registry = load_sequence_registry(sequence_registry_path)

    host_genome = _resolve_sequence(
        registry,
        kind="host_genome",
        expected_sequence_sha256=resources.host_state["host_genome_sha256"],
    )

    sequence = resources.phage_profile["sequence"]
    if sequence is None:
        raise ValueError(
            "PhageProfile has no sequence reference; sequence resolution cannot continue."
        )
    if not isinstance(sequence, dict):
        raise ValueError("PhageProfile sequence reference must be an object.")
    expected_phage_sha256 = sequence.get("sha256")
    if not _is_sha256(expected_phage_sha256):
        raise ValueError(
            "PhageProfile sequence reference must include a normalized SHA-256; "
            "accession-only references cannot be resolved exactly."
        )

    phage_sequence = _resolve_sequence(
        registry,
        kind="phage_sequence",
        expected_sequence_sha256=expected_phage_sha256,
        expected_accession=sequence.get("accession"),
        expected_length=sequence.get("length"),
    )

    return ResolvedModelSequences(
        resources=resources,
        sequence_registry_id=registry.registry_id,
        sequence_registry_file_sha256=registry.file_sha256,
        sequence_registry_root=registry.root,
        host_genome=host_genome,
        phage_sequence=phage_sequence,
    )
