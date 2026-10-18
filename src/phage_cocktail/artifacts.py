"""Shared integrity helpers for versioned JSON artifacts and local files."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def read_json_object(path: str | Path, *, label: str = "JSON artifact") -> dict:
    """Read a JSON object while rejecting non-standard NaN/Infinity constants."""
    source = Path(path)

    def reject_constant(value: str):
        raise ValueError(f"non-standard JSON constant {value!r} is not allowed")

    try:
        data = json.loads(
            source.read_text(encoding="utf-8"),
            parse_constant=reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
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
