"""Shared integrity helpers for versioned JSON artifacts and local files."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

PACKAGE_DIR = Path(__file__).resolve().parent


def read_json_object(path: str | Path, *, label: str = "JSON artifact") -> dict:
    """Read a JSON object while rejecting non-standard NaN/Infinity constants."""
    source = Path(path)

    def reject_constant(value: str):
        raise ValueError(f"non-standard JSON constant {value!r} is not allowed")

    try:
        data = json.loads(source.read_text(encoding="utf-8"), parse_constant=reject_constant)
    except (OSError, UnicodeError, ValueError) as exc:
        raise ValueError(f"cannot read valid JSON for {label} from {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{label} at {source} must contain a JSON object.")
    return data


def canonical_json(document: dict) -> str:
    """Serialize JSON deterministically using the model/lab exchange convention."""
    return json.dumps(
        document,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def canonical_content_sha256(document: dict) -> str:
    """Hash canonical record content, excluding an optional content_sha256 field."""
    if not isinstance(document, dict):
        raise ValueError("content to hash must be a JSON object.")
    content = dict(document)
    content.pop("content_sha256", None)
    try:
        canonical = canonical_json(content).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"document is not canonical JSON data: {exc}") from exc
    return hashlib.sha256(canonical).hexdigest()


def sha256_file(path: str | Path) -> str:
    """Hash the exact bytes of a local file without loading the whole file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def code_provenance(repo: str | Path, script: str | Path) -> dict:
    """Record the Git revision and hashes of the entry script and every package module.

    Hashing the whole package avoids hand-maintained file lists that drift from
    what a run actually imports.
    """
    repo = Path(repo).resolve()
    paths = [Path(script).resolve(), *sorted(PACKAGE_DIR.glob("*.py"))]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    return {
        "code_base_revision": revision,
        "code_sha256": {path.relative_to(repo).as_posix(): sha256_file(path) for path in paths},
    }


@contextmanager
def staged_directory(destination: str | Path):
    """Yield a sibling temporary directory that becomes ``destination`` on success."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        yield stage
        stage.rename(destination)
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def atomic_write_json(document: dict, path: str | Path) -> None:
    destination = Path(path)
    temporary_name = None
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            json.dump(document, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, destination)
    except (OSError, TypeError, ValueError) as exc:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except OSError:
                pass
        raise ValueError(f"cannot write JSON artifact to {destination}: {exc}") from exc
