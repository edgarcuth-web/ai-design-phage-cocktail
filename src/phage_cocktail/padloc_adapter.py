"""Run pinned PADLOC against a registry-verified host genome.

The adapter invokes the external PADLOC executable; it never fabricates feature
calls. A successful zero-hit run is represented explicitly from PADLOC's own
completion message. Tool/database versions are pinned to v2.0.0 for this
adapter version, and the installed database tree is content-addressed.
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tempfile

from jsonschema import Draft202012Validator, FormatChecker

from .artifacts import (
    canonical_content_sha256,
    canonical_json,
    read_json_object,
    sha256_file,
)
from .feature_artifact import build_feature_artifact, write_feature_artifact
from .resource_registry import load_host_state_record
from .sequence_registry import resolve_registered_host_genome

PADLOC_VERSION = "v2.0.0"
PADLOC_DB_VERSION = "v2.0.0"
PADLOC_DB_RELEASE_URI = "https://github.com/padlocbio/padloc-db/releases/tag/v2.0.0"
PADLOC_DB_HASH_METHOD = "sha256-canonical-json-of-sorted-relative-file-digests-v1"
PADLOC_RESULT_SCHEMA_PATH = (
    Path(__file__).resolve().parent
    / "schemas"
    / "padloc-host-defense-result-v1.schema.json"
)
DEFAULT_TIMEOUT_SECONDS = 7200
_REQUIRED_DATABASE_FILES = (
    "RELEASE.md",
    "hmm/padlocdb.hmm",
    "cm/padlocdb.cm",
    "hmm_meta.txt",
    "sys_meta.txt",
)
_REQUIRED_CSV_COLUMNS = ("system.number", "seqid", "system", "target.name")
_LIMITATIONS = [
    "PADLOC searches known defense systems represented in PADLOC-DB v2.0.0; no hits do not prove that the genome lacks anti-phage defense systems.",
    "The --fna workflow predicts protein-coding genes with Prodigal. This run did not supply separate CRISPR-array or ncRNA annotations, so those feature classes may be incomplete.",
    "The result is computational annotation, not experimental validation. Review native PADLOC outputs and source evidence before admitting calls into HostProfile.",
    "PADLOC and PADLOC-DB are pinned; transitive dependency versions are recorded in the dependency-check log and environment metadata rather than all being independently pinned by this adapter.",
]


@lru_cache(maxsize=1)
def _result_validator() -> Draft202012Validator:
    schema = read_json_object(
        PADLOC_RESULT_SCHEMA_PATH, label="PADLOC host-defense result schema"
    )
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, format_checker=FormatChecker())


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def _capture(command: list[str], *, timeout_seconds: int = 60) -> str:
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout_seconds,
            check=False,
        )
    except FileNotFoundError as exc:
        raise ValueError(f"required executable was not found: {command[0]!r}.") from exc
    except subprocess.TimeoutExpired as exc:
        raise ValueError(
            f"command timed out after {timeout_seconds} seconds: {command[0]!r}."
        ) from exc
    except OSError as exc:
        raise ValueError(f"cannot execute {command[0]!r}: {exc}") from exc
    output = (completed.stdout or "") + (completed.stderr or "")
    if completed.returncode != 0:
        tail = "\n".join(output.splitlines()[-12:])
        raise ValueError(
            f"command failed with exit code {completed.returncode}: "
            f"{command[0]!r}. Output tail:\n{tail}"
        )
    return output


def _resolve_executable(executable: str | Path) -> Path:
    candidate_text = str(executable)
    located = shutil.which(candidate_text)
    if located is None:
        candidate = Path(candidate_text).expanduser()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            located = str(candidate)
    if located is None:
        raise ValueError(
            f"PADLOC executable {candidate_text!r} was not found on PATH; "
            "install the pinned environment described in docs/padloc-adapter.md."
        )
    resolved = Path(located).resolve(strict=True)
    if not resolved.is_file() or not os.access(resolved, os.X_OK):
        raise ValueError(f"PADLOC executable is not an executable file: {resolved}")
    return resolved


def _parse_version(output: str, *, pattern: str, label: str) -> str:
    match = re.search(pattern, output, flags=re.IGNORECASE | re.MULTILINE)
    if match is None:
        raise ValueError(f"could not parse {label} from command output: {output.strip()!r}")
    return match.group(1)


def _database_snapshot(data_dir: str | Path) -> dict:
    root = Path(data_dir).expanduser()
    try:
        root = root.resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"PADLOC data directory not found: {data_dir}") from exc
    if not root.is_dir():
        raise ValueError(f"PADLOC data path is not a directory: {root}")

    for relative in _REQUIRED_DATABASE_FILES:
        path = root / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(
                f"PADLOC-DB v2.0.0 is incomplete: required file is missing or empty: {relative}"
            )
    release_text = (root / "RELEASE.md").read_text(encoding="utf-8")
    if "# PADLOC-DB release v2.0.0" not in release_text.splitlines():
        raise ValueError(
            "PADLOC-DB release marker is not exactly v2.0.0; refusing an unpinned database."
        )
    systems_dir = root / "sys"
    if not systems_dir.is_dir() or not any(systems_dir.glob("*.yaml")):
        raise ValueError("PADLOC-DB v2.0.0 has no system model YAML files.")

    entries = []
    total_bytes = 0
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        directory_path = Path(directory)
        for name in sorted(dirnames):
            child = directory_path / name
            if child.is_symlink():
                raise ValueError(f"PADLOC-DB snapshot must not contain symlink directories: {child}")
        for name in sorted(filenames):
            path = directory_path / name
            if path.is_symlink():
                raise ValueError(f"PADLOC-DB snapshot must not contain symlink files: {path}")
            if not path.is_file():
                continue
            digest = sha256_file(path)
            size = path.stat().st_size
            entries.append(
                {
                    "path": path.relative_to(root).as_posix(),
                    "sha256": digest,
                    "size_bytes": size,
                }
            )
            total_bytes += size
    entries.sort(key=lambda item: item["path"])
    if not entries:
        raise ValueError("PADLOC-DB snapshot contains no regular files.")
    tree_sha256 = hashlib.sha256(canonical_json(entries).encode("utf-8")).hexdigest()
    return {
        "root": root,
        "version": PADLOC_DB_VERSION,
        "tree_sha256": tree_sha256,
        "hash_method": PADLOC_DB_HASH_METHOD,
        "file_count": len(entries),
        "total_bytes": total_bytes,
    }


def _prepare_output_directory(
    path: str | Path, *, database_root: Path, sequence_registry_root: Path
) -> Path:
    requested = Path(path).expanduser()
    if requested.is_symlink():
        raise ValueError("PADLOC output directory must not itself be a symlink.")
    root = requested.resolve(strict=False)
    if root == database_root or root.is_relative_to(database_root) or database_root.is_relative_to(root):
        raise ValueError("PADLOC output directory and database directory must not overlap.")
    if root.is_relative_to(sequence_registry_root):
        raise ValueError(
            "PADLOC output directory must not be inside the sequence registry directory."
        )
    if root.exists():
        if not root.is_dir():
            raise ValueError(f"PADLOC output path is not a directory: {root}")
        if any(root.iterdir()):
            raise ValueError(
                f"PADLOC output directory must be empty to prevent mixing runs: {root}"
            )
    else:
        root.mkdir(parents=True, exist_ok=True)
    return root


def _parse_native_csv(path: Path) -> list[dict]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames:
                raise ValueError("PADLOC output CSV has no header.")
            missing = sorted(set(_REQUIRED_CSV_COLUMNS) - set(reader.fieldnames))
            if missing:
                raise ValueError(
                    "PADLOC output CSV is missing required columns: "
                    + ", ".join(missing)
                )
            records = []
            for row_number, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(
                        f"PADLOC output CSV row {row_number} has more values than its header."
                    )
                if any(row.get(field) is None or not row[field].strip() for field in _REQUIRED_CSV_COLUMNS):
                    raise ValueError(
                        f"PADLOC output CSV row {row_number} has a blank required field."
                    )
                records.append(dict(row))
            return records
    except OSError as exc:
        raise ValueError(f"cannot read PADLOC output CSV {path}: {exc}") from exc


def validate_padloc_summary(document: dict) -> dict:
    """Validate summary schema, canonical digest, and hit-count semantics."""
    if not isinstance(document, dict):
        raise ValueError("PADLOC result summary must be a JSON object.")
    errors = sorted(
        _result_validator().iter_errors(document),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        summaries = []
        for error in errors[:5]:
            location = ".".join(str(part) for part in error.absolute_path) or "$"
            summaries.append(f"{location}: {error.message}")
        extra = f" (and {len(errors) - 5} more)" if len(errors) > 5 else ""
        raise ValueError("invalid PADLOC result summary: " + "; ".join(summaries) + extra)

    expected_hash = canonical_content_sha256(document)
    if document["content_sha256"] != expected_hash:
        raise ValueError(
            "PADLOC result summary content_sha256 mismatch: "
            f"expected {expected_hash}, got {document['content_sha256']!r}."
        )
    records = document["records"]
    if document["gene_hit_count"] != len(records):
        raise ValueError("PADLOC result summary gene_hit_count does not match records.")
    systems = {
        (row["seqid"], row["system.number"], row["system"])
        for row in records
    }
    system_names = sorted({row["system"] for row in records})
    if document["systems_detected"] != len(systems):
        raise ValueError("PADLOC result summary systems_detected does not match records.")
    if document["system_names"] != system_names:
        raise ValueError("PADLOC result summary system_names does not match records.")
    expected_outcome = "hits_detected" if records else "no_hits_detected"
    if document["outcome"] != expected_outcome:
        raise ValueError("PADLOC result summary outcome does not match records.")
    if records and document["native_output_path"] is None:
        raise ValueError("PADLOC result summary with hits must retain the native CSV path.")
    return document


def load_padloc_summary(path: str | Path) -> dict:
    """Load a PADLOC summary and validate its schema and canonical digest."""
    return validate_padloc_summary(
        read_json_object(path, label="PADLOC host-defense result")
    )


def _native_output_specs(native_dir: Path, artifact_root: Path) -> list[dict]:
    formats = {
        ".csv": ("csv", "text/csv"),
        ".gff": ("gff3", "text/gff3"),
        ".faa": ("fasta", "text/x-fasta"),
        ".domtblout": ("hmmer-domtblout", "text/plain"),
    }
    specs = []
    for path in sorted(native_dir.iterdir(), key=lambda item: item.name):
        if path.is_symlink():
            raise ValueError(f"PADLOC created an unexpected symlink output: {path.name}")
        if not path.is_file():
            raise ValueError(f"PADLOC created an unexpected non-file output: {path.name}")
        file_format, media_type = formats.get(
            path.suffix, (path.suffix.lstrip(".") or "binary", "application/octet-stream")
        )
        specs.append(
            {
                "role": "raw",
                "path": path.relative_to(artifact_root).as_posix(),
                "format": file_format,
                "media_type": media_type,
            }
        )
    return specs


def run_padloc_host_defense(
    *,
    host_state_path: str | Path,
    sequence_registry_path: str | Path,
    padloc_data_dir: str | Path,
    output_dir: str | Path,
    padloc_executable: str | Path = "padloc",
    cpu: int = 1,
    timeout_seconds: int = DEFAULT_TIMEOUT_SECONDS,
) -> dict:
    """Run actual PADLOC v2.0.0 and emit a verified FeatureArtifact plus summary.

    PADLOC and PADLOC-DB are version-pinned. The input must be a valid HostState
    v1 record whose exact genome fingerprint resolves through the local sequence
    registry. The caller must install PADLOC-DB v2.0.0 in a dedicated data
    directory before invoking this function; this function never downloads or
    mutates the database.
    """
    if type(cpu) is not int or cpu < 1:
        raise ValueError("cpu must be a positive integer.")
    if type(timeout_seconds) is not int or timeout_seconds < 1:
        raise ValueError("timeout_seconds must be a positive integer.")

    host_state, host_state_file_sha256 = load_host_state_record(host_state_path)
    resolved_genome = resolve_registered_host_genome(
        sequence_registry_path, host_state["host_genome_sha256"]
    )
    asset = resolved_genome.asset

    executable = _resolve_executable(padloc_executable)
    version_output = _capture([str(executable), "--version"])
    observed_version = _parse_version(
        version_output,
        pattern=r"^\s*padloc\s+(v?\d+\.\d+\.\d+)\s*$",
        label="PADLOC version",
    )
    if observed_version != PADLOC_VERSION:
        raise ValueError(
            f"PADLOC version mismatch: expected {PADLOC_VERSION}, got {observed_version}."
        )

    database_root_arg = Path(padloc_data_dir).expanduser()
    try:
        database_root = database_root_arg.resolve(strict=True)
    except OSError as exc:
        raise ValueError(f"PADLOC data directory not found: {padloc_data_dir}") from exc
    db_output = _capture(
        [str(executable), "--data", str(database_root), "--db-version"]
    )
    observed_db_version = _parse_version(
        db_output,
        pattern=r"^\s*padloc-db\s+(v?\d+\.\d+\.\d+)\s*$",
        label="PADLOC-DB version",
    )
    if observed_db_version != PADLOC_DB_VERSION:
        raise ValueError(
            f"PADLOC-DB version mismatch: expected {PADLOC_DB_VERSION}, "
            f"got {observed_db_version}."
        )
    database = _database_snapshot(database_root)

    artifact_root = _prepare_output_directory(
        output_dir,
        database_root=database_root,
        sequence_registry_root=resolved_genome.sequence_registry_root,
    )
    logs_dir = artifact_root / "logs"
    native_dir = artifact_root / "outputs" / "padloc-native"
    logs_dir.mkdir(parents=True, exist_ok=True)
    native_dir.mkdir(parents=True, exist_ok=True)
    _write_text(logs_dir / "padloc-version.txt", version_output)
    _write_text(logs_dir / "padloc-db-version.txt", db_output)

    check_output = _capture(
        [str(executable), "--check-deps"], timeout_seconds=120
    )
    _write_text(logs_dir / "padloc-check-deps.txt", check_output)
    if re.search(r"\bWARNING\b", check_output, flags=re.IGNORECASE):
        raise ValueError(
            "PADLOC dependency check reported missing dependencies; "
            f"inspect {logs_dir / 'padloc-check-deps.txt'}."
        )

    temporary_root = Path("/tmp") if Path("/tmp").is_dir() else Path(tempfile.gettempdir())
    if any(character.isspace() for character in str(temporary_root.resolve())):
        raise ValueError(
            "PADLOC requires a temporary input path without whitespace because "
            "the upstream wrapper expands the FASTA basename without quoting it."
        )

    started_at = _timestamp()
    stdout = ""
    stderr = ""
    with tempfile.TemporaryDirectory(
        prefix="padloc-host-", dir=str(temporary_root.resolve())
    ) as temporary_directory:
        input_copy = Path(temporary_directory) / "host.fna"
        try:
            shutil.copyfile(asset.path, input_copy)
        except OSError as exc:
            raise ValueError(f"cannot stage verified host FASTA for PADLOC: {exc}") from exc
        if sha256_file(input_copy) != asset.file_sha256:
            raise ValueError("staged host FASTA does not match the registered input file digest.")

        command = [
            str(executable),
            "--data",
            str(database_root),
            "--fna",
            str(input_copy),
            "--outdir",
            str(native_dir),
            "--cpu",
            str(cpu),
        ]
        try:
            run = subprocess.run(
                command,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout_seconds,
                check=False,
            )
            stdout = run.stdout or ""
            stderr = run.stderr or ""
        except subprocess.TimeoutExpired as exc:
            stdout = exc.stdout.decode("utf-8", "replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            stderr = exc.stderr.decode("utf-8", "replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
            _write_text(logs_dir / "padloc.stdout.txt", stdout)
            _write_text(logs_dir / "padloc.stderr.txt", stderr)
            raise ValueError(
                f"PADLOC exceeded the {timeout_seconds}-second timeout; "
                f"inspect {logs_dir} for partial logs."
            ) from exc
        except OSError as exc:
            raise ValueError(f"cannot execute PADLOC: {exc}") from exc

    finished_at = _timestamp()
    _write_text(logs_dir / "padloc.stdout.txt", stdout)
    _write_text(logs_dir / "padloc.stderr.txt", stderr)
    if run.returncode != 0:
        tail = "\n".join((stdout + stderr).splitlines()[-15:])
        raise ValueError(
            f"PADLOC exited with code {run.returncode}. Logs retained in {logs_dir}. "
            f"Output tail:\n{tail}"
        )

    combined_output = stdout + "\n" + stderr
    native_csv = native_dir / "host_padloc.csv"
    no_hits_message = re.search(
        r"Complete:\s*Nothing found for\s+host\b",
        combined_output,
        flags=re.IGNORECASE,
    )
    if native_csv.is_file():
        records = _parse_native_csv(native_csv)
        if not records and no_hits_message is None:
            raise ValueError(
                "PADLOC produced a zero-row CSV without its explicit no-hit completion message; "
                f"logs retained in {logs_dir}."
            )
    elif no_hits_message is not None:
        records = []
    else:
        raise ValueError(
            "PADLOC returned success but produced neither its native CSV nor an explicit "
            f"no-hit completion message; logs retained in {logs_dir}."
        )

    systems = {
        (row["seqid"], row["system.number"], row["system"])
        for row in records
    }
    system_names = sorted({row["system"] for row in records})
    run_stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    artifact_id = f"padloc-host-defense:{host_state['record_id']}:{run_stamp}"
    native_output_relative = (
        native_csv.relative_to(artifact_root).as_posix() if native_csv.is_file() else None
    )
    summary = {
        "schema_version": 1,
        "artifact_id": artifact_id,
        "created_at": _timestamp(),
        "subject": {
            "id": host_state["record_id"],
            "record_sha256": host_state["content_sha256"],
            "record_file_sha256": host_state_file_sha256,
            "sequence_sha256": asset.sequence_sha256,
            "input_file_sha256": asset.file_sha256,
            "length_bases": asset.bases,
            "contigs": asset.contigs,
            "sequence_registry_id": resolved_genome.sequence_registry_id,
            "sequence_registry_file_sha256": resolved_genome.sequence_registry_file_sha256,
        },
        "analyzer": {
            "name": "PADLOC",
            "version": observed_version,
            "executable_sha256": sha256_file(executable),
            "database": {
                "name": "PADLOC-DB",
                "version": observed_db_version,
                "tree_sha256": database["tree_sha256"],
                "hash_method": database["hash_method"],
                "file_count": database["file_count"],
                "total_bytes": database["total_bytes"],
                "source_uri": PADLOC_DB_RELEASE_URI,
            },
        },
        "outcome": "hits_detected" if records else "no_hits_detected",
        "systems_detected": len(systems),
        "gene_hit_count": len(records),
        "system_names": system_names,
        "records": records,
        "native_output_path": native_output_relative,
        "limitations": list(_LIMITATIONS),
    }
    summary["content_sha256"] = canonical_content_sha256(summary)
    validate_padloc_summary(summary)
    summary_path = artifact_root / "outputs" / "host-defense-summary.json"
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    log_specs = [
        {
            "role": "log",
            "path": path.relative_to(artifact_root).as_posix(),
            "format": "txt",
            "media_type": "text/plain",
        }
        for path in sorted(logs_dir.iterdir(), key=lambda item: item.name)
        if path.is_file()
    ]
    output_specs = [
        {
            "role": "primary",
            "path": summary_path.relative_to(artifact_root).as_posix(),
            "format": "json",
            "media_type": "application/json",
            "schema_sha256": sha256_file(PADLOC_RESULT_SCHEMA_PATH),
        },
        *_native_output_specs(native_dir, artifact_root),
        *log_specs,
    ]

    artifact = build_feature_artifact(
        manifest_path=artifact_root / "feature-artifact.v1.json",
        artifact_id=artifact_id,
        feature_set="host_defense_systems",
        subject={
            "kind": "host_genome",
            "id": host_state["record_id"],
            "sequence_sha256": asset.sequence_sha256,
            "input_file_sha256": asset.file_sha256,
            "length_bases": asset.bases,
            "contigs": asset.contigs,
            "sequence_registry_id": resolved_genome.sequence_registry_id,
            "sequence_registry_file_sha256": resolved_genome.sequence_registry_file_sha256,
            "source_record_sha256": host_state["content_sha256"],
            "source_record_file_sha256": host_state_file_sha256,
        },
        tool={
            "name": "PADLOC",
            "version": observed_version,
            "executable_sha256": sha256_file(executable),
        },
        parameters={
            "input_mode": "--fna",
            "cpu": cpu,
            "timeout_seconds": timeout_seconds,
            "padloc_version": observed_version,
            "padloc_db_version": observed_db_version,
            "padloc_db_tree_sha256": database["tree_sha256"],
            "padloc_db_hash_method": database["hash_method"],
            "padloc_db_file_count": database["file_count"],
            "padloc_db_total_bytes": database["total_bytes"],
            "gene_calling": "PADLOC default Prodigal workflow",
            "crisper_array_input_supplied": False,
            "ncrna_input_supplied": False,
        },
        reference_data=[
            {
                "name": "PADLOC-DB",
                "version": observed_db_version,
                "sha256": database["tree_sha256"],
                "source_uri": PADLOC_DB_RELEASE_URI,
            }
        ],
        outputs=output_specs,
        started_at=started_at,
        finished_at=finished_at,
        created_at=_timestamp(),
        status="complete",
        environment={
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "padloc_version_output_sha256": sha256_file(logs_dir / "padloc-version.txt"),
            "padloc_db_version_output_sha256": sha256_file(logs_dir / "padloc-db-version.txt"),
            "dependency_check_output_sha256": sha256_file(logs_dir / "padloc-check-deps.txt"),
            "dependency_check_output": check_output.strip(),
            "requested_cpu_count": cpu,
            "database_snapshot_file_count": database["file_count"],
            "database_snapshot_total_bytes": database["total_bytes"],
            "command_template": [
                "padloc",
                "--data",
                "<PADLOC_DB_ROOT>",
                "--fna",
                "<verified-host-genome>",
                "--outdir",
                "<FeatureArtifact>/outputs/padloc-native",
                "--cpu",
                str(cpu),
            ],
        },
        limitations=list(_LIMITATIONS),
    )
    manifest_path = artifact_root / "feature-artifact.v1.json"
    write_feature_artifact(artifact, manifest_path, verify_outputs=True)
    return {
        "manifest_path": manifest_path,
        "artifact": artifact,
        "summary_path": summary_path,
        "summary": summary,
        "dependency_check": check_output,
        "database": database,
    }
