"""Shared integrity helpers for versioned JSON artifacts and local files."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path


def read_json_object_with_digest(
    path: str | Path,
    *,
    label: str = "JSON artifact",
) -> tuple[dict, str]:
    """Read a JSON object and hash the exact same bytes to avoid read/hash races."""
    source = Path(path)

    def reject_constant(value: str):
        raise ValueError(f"non-standard JSON constant {value!r} is not allowed")

    try:
        raw = source.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        data = json.loads(
            raw.decode("utf-8"),
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"cannot read valid JSON for {label} from {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{label} at {source} must contain a JSON object.")
    return data, digest


def read_json_object(path: str | Path, *, label: str = "JSON artifact") -> dict:
    """Read a JSON object while rejecting non-standard NaN/Infinity constants."""
    document, _ = read_json_object_with_digest(path, label=label)
    return document

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
