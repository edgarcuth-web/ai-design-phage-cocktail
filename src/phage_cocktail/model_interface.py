"""Strict file-based I/O for the external ModelRequest/ModelResult v1 boundary.

This module validates exchange artifacts and records explicit unresolved results.
It does not resolve biological features, train models, or run inference.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from functools import lru_cache
from collections.abc import Mapping

from jsonschema import Draft202012Validator, FormatChecker

SCHEMA_VERSION = 1
PREDICTION_TARGETS = (
    "adsorption",
    "productive_infection",
    "resistance_emergence",
)
SCHEMA_DIR = Path(__file__).resolve().parents[2] / "schemas"
UNKNOWN_ESTIMATE = {
    "status": "not_estimated",
    "probability": None,
    "uncertainty_lower": None,
    "uncertainty_upper": None,
    "uncertainty_method": None,
}


def _reject_json_constant(value: str):
    raise ValueError(f"non-standard JSON constant {value!r} is not allowed")


def _read_json(path: str | Path) -> dict:
    source = Path(path)
    try:
        data = json.loads(
            source.read_text(encoding="utf-8"),
            parse_constant=_reject_json_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"cannot read valid JSON from {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{source} must contain a JSON object.")
    return data


def _canonical_json(data: dict) -> str:
    return json.dumps(
        data,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _content_hash(document: dict) -> str:
    content = dict(document)
    content.pop("content_sha256", None)
    try:
        canonical = _canonical_json(content).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"document is not canonical JSON data: {exc}") from exc
    return hashlib.sha256(canonical).hexdigest()


def _verify_content_hash(document: dict, kind: str) -> None:
    supplied = document.get("content_sha256")
    expected = _content_hash(document)
    if supplied != expected:
        raise ValueError(
            f"{kind} content_sha256 mismatch: expected {expected}, got {supplied!r}."
        )


def _parse_timestamp(value: str, field_name: str) -> datetime:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.") from None
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{field_name} must include a timezone.")
    return parsed


@lru_cache(maxsize=2)
def _validator(schema_name: str) -> Draft202012Validator:
    if schema_name not in {"model-request-v1.schema.json", "model-result-v1.schema.json"}:
        raise ValueError(f"unsupported schema name: {schema_name}")
    schema = _read_json(SCHEMA_DIR / schema_name)
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, format_checker=FormatChecker())


def _validate_schema(document: dict, schema_name: str, kind: str) -> None:
    if not isinstance(document, dict):
        raise ValueError(f"{kind} must be a JSON object.")
    errors = sorted(
        _validator(schema_name).iter_errors(document),
        key=lambda error: tuple(str(part) for part in error.absolute_path),
    )
    if errors:
        summaries = []
        for error in errors[:5]:
            location = ".".join(str(part) for part in error.absolute_path) or "$"
            summaries.append(f"{location}: {error.message}")
        extra = f" (and {len(errors) - 5} more)" if len(errors) > 5 else ""
        raise ValueError(f"invalid {kind}: " + "; ".join(summaries) + extra)


def validate_model_request(request: dict) -> dict:
    """Validate ModelRequest schema, content hash, derived ID and timestamp order."""
    _validate_schema(request, "model-request-v1.schema.json", "ModelRequest")
    _verify_content_hash(request, "ModelRequest")
    if request["record_id"] != request["request_id"]:
        raise ValueError("ModelRequest record_id must equal request_id.")
    if _parse_timestamp(request["created_at"], "created_at") < _parse_timestamp(
        request["host_state"]["as_of"], "host_state.as_of"
    ):
        raise ValueError("ModelRequest created_at cannot predate HostState as_of.")
    return request


def load_model_request(path: str | Path) -> dict:
    """Read and fully validate a ModelRequest JSON file."""
    return validate_model_request(_read_json(path))


def validate_model_result(result: dict, request: dict | None = None) -> dict:
    """Validate ModelResult schema/hash and, optionally, exact request lineage."""
    _validate_schema(result, "model-result-v1.schema.json", "ModelResult")
    _verify_content_hash(result, "ModelResult")
    if result["record_id"] != result["result_id"]:
        raise ValueError("ModelResult record_id must equal result_id.")

    predicted_at = _parse_timestamp(result["predicted_at"], "predicted_at")
    host_as_of = _parse_timestamp(result["host_state"]["as_of"], "host_state.as_of")
    if predicted_at < host_as_of:
        raise ValueError("ModelResult predicted_at cannot predate HostState as_of.")

    for target in PREDICTION_TARGETS:
        estimate = result["probabilities"][target]
        if estimate["status"] == "estimated":
            if not (
                estimate["uncertainty_lower"]
                <= estimate["probability"]
                <= estimate["uncertainty_upper"]
            ):
                raise ValueError(
                    f"{target} uncertainty interval must contain its probability."
                )
            if result["producer"]["artifact_sha256"] is None:
                raise ValueError(
                    "estimated output requires producer.artifact_sha256."
                )

    if request is not None:
        request = validate_model_request(request)
        if result["request_id"] != request["request_id"]:
            raise ValueError("ModelResult request_id does not match ModelRequest.")
        if result["request_sha256"] != request["content_sha256"]:
            raise ValueError("ModelResult request_sha256 does not match exact ModelRequest.")
        if result["host_state"] != request["host_state"]:
            raise ValueError("ModelResult HostState does not match ModelRequest.")
        if result["phage_profile_id"] != request["phage_profile_id"]:
            raise ValueError("ModelResult PhageProfile ID does not match ModelRequest.")
        if result["phage_profile_sha256"] != request["phage_profile_sha256"]:
            raise ValueError("ModelResult PhageProfile fingerprint does not match ModelRequest.")
        if predicted_at < _parse_timestamp(request["created_at"], "request.created_at"):
            raise ValueError("ModelResult predicted_at cannot predate ModelRequest creation.")
        requested = set(request["requested_targets"])
        for target in PREDICTION_TARGETS:
            if target not in requested and result["probabilities"][target]["status"] == "estimated":
                raise ValueError(f"ModelResult estimated unrequested target {target!r}.")
    return result


def load_model_result(path: str | Path, request: dict | None = None) -> dict:
    """Read and validate a ModelResult JSON file, optionally against its request."""
    return validate_model_result(_read_json(path), request=request)


def build_model_result(
    request: dict,
    *,
    result_id: str,
    probabilities: Mapping[str, dict],
    producer: dict,
    mechanistic_explanation: str,
    limitations: tuple[str, ...] | list[str],
    predicted_at: str | None = None,
) -> dict:
    """Build a content-addressed result from explicit model outputs.

    Missing target entries become not_estimated; unknown target names are rejected.
    This function never converts missing observations into negative labels.
    """
    request = validate_model_request(request)
    if not isinstance(probabilities, Mapping):
        raise ValueError("probabilities must be a mapping of target names to estimates.")
    unsupported = set(probabilities) - set(PREDICTION_TARGETS)
    if unsupported:
        raise ValueError(
            "unsupported prediction target(s): " + ", ".join(sorted(map(str, unsupported)))
        )
    if not isinstance(producer, dict):
        raise ValueError("producer must be an object.")
    if not isinstance(limitations, (tuple, list)):
        raise ValueError("limitations must be a list or tuple of strings.")
    if predicted_at is None:
        predicted_at = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")

    result = {
        "schema_version": SCHEMA_VERSION,
        "record_id": result_id,
        "result_id": result_id,
        "request_id": request["request_id"],
        "request_sha256": request["content_sha256"],
        "host_state": request["host_state"],
        "phage_profile_id": request["phage_profile_id"],
        "phage_profile_sha256": request["phage_profile_sha256"],
        "probabilities": {
            target: dict(probabilities.get(target, UNKNOWN_ESTIMATE))
            for target in PREDICTION_TARGETS
        },
        "mechanistic_explanation": mechanistic_explanation,
        "limitations": list(limitations),
        "producer": dict(producer),
        "predicted_at": predicted_at,
    }
    result["content_sha256"] = _content_hash(result)
    return validate_model_result(result, request=request)


def create_unresolved_result(
    request: dict,
    *,
    result_id: str | None = None,
    predicted_at: str | None = None,
) -> dict:
    """Create a valid, explicitly non-predictive response when inference is unavailable."""
    request = validate_model_request(request)
    if result_id is None:
        result_id = f"{request['request_id']}:unresolved"
    return build_model_result(
        request,
        result_id=result_id,
        probabilities={},
        producer={
            "name": "unresolved-no-inference",
            "version": "1",
            "parameters": {"reason": "inference_not_available"},
            "artifact_sha256": None,
        },
        mechanistic_explanation=(
            "No inference was run. This result records request lineage only and "
            "contains no biological prediction."
        ),
        limitations=(
            "No inference implementation was run; all targets remain not_estimated.",
            "The request contains identifiers and fingerprints, not resolved biological features.",
            "This software response is not biological evidence and must not be used as a training label.",
        ),
        predicted_at=predicted_at,
    )


def _atomic_write_json(document: dict, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name = None
    try:
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


def write_model_result(
    result: dict,
    path: str | Path,
    request: dict | None = None,
) -> None:
    """Validate before atomically writing a ModelResult JSON file."""
    validate_model_result(result, request=request)
    _atomic_write_json(result, path)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phage-model-interface",
        description="Validate v1 exchange files; this CLI does not run biological inference.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    validate_request = commands.add_parser(
        "validate-request", help="validate a ModelRequest JSON file"
    )
    validate_request.add_argument("request", type=Path)

    validate_result = commands.add_parser(
        "validate-result", help="validate a ModelResult JSON file"
    )
    validate_result.add_argument("result", type=Path)
    validate_result.add_argument(
        "--request",
        type=Path,
        help="also verify exact request/result lineage and requested-target scope",
    )

    unresolved = commands.add_parser(
        "create-unresolved-result",
        help="write a not_estimated response without running inference",
    )
    unresolved.add_argument("request", type=Path)
    unresolved.add_argument("output", type=Path)
    unresolved.add_argument("--result-id")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate-request":
            request = load_model_request(args.request)
            print(
                f"valid ModelRequest v{request['schema_version']}: "
                f"{request['request_id']} sha256={request['content_sha256']}"
            )
            return 0

        if args.command == "validate-result":
            request = load_model_request(args.request) if args.request else None
            result = load_model_result(args.result, request=request)
            suffix = " with exact request lineage" if request is not None else ""
            print(
                f"valid ModelResult v{result['schema_version']}{suffix}: "
                f"{result['result_id']} sha256={result['content_sha256']}"
            )
            return 0

        if args.command == "create-unresolved-result":
            if args.request.resolve() == args.output.resolve():
                raise ValueError("output path must differ from request path.")
            request = load_model_request(args.request)
            result = create_unresolved_result(request, result_id=args.result_id)
            write_model_result(result, args.output, request=request)
            print(
                f"wrote unresolved ModelResult: {args.output} "
                f"sha256={result['content_sha256']}; no inference was run"
            )
            return 0

        parser.error(f"unknown command: {args.command}")
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
