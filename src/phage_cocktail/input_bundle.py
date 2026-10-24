"""Consume portable SDL exports without importing laboratory implementation.

This module verifies wire integrity. SDL is responsible for admitting scientific
snapshots and resolving sequence identity before export. Models transform these
frozen inputs; they never fetch or annotate missing inputs implicitly.
"""

from datetime import datetime, timezone
from functools import lru_cache
import hashlib
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker

from .artifacts import canonical_content_sha256, canonical_json, read_json_object

SCHEMA_DIR = Path(__file__).resolve().parent / "schemas"
SCHEMAS = {"model-input-bundle-v1.schema.json", "pair-prediction-v1.schema.json"}


def _time(value, field_name="timestamp"):
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        raise ValueError(f"{field_name} must be an ISO-8601 timestamp.") from None
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{field_name} must include a timezone.")
    return parsed


@lru_cache(maxsize=len(SCHEMAS))
def _validator(schema_name: str) -> Draft202012Validator:
    if schema_name not in SCHEMAS:
        raise ValueError(f"unsupported schema name: {schema_name}")
    schema = read_json_object(SCHEMA_DIR / schema_name)
    Draft202012Validator.check_schema(schema)
    return Draft202012Validator(schema, format_checker=FormatChecker())


def _validate_schema(document, schema_name: str, kind: str) -> None:
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


def _verify_content_hash(document: dict, kind: str) -> None:
    supplied = document.get("content_sha256")
    expected = canonical_content_sha256(document)
    if supplied != expected:
        raise ValueError(f"{kind} content_sha256 mismatch: expected {expected}, got {supplied!r}.")


def _canonical_sha256(value) -> str:
    return hashlib.sha256(canonical_json(value).encode()).hexdigest()


def _file(root, relative, expected_hash, expected_size):
    if (not relative or "\\" in relative or ":" in relative
            or any(p in {"", ".", ".."} for p in relative.split("/"))):
        raise ValueError("resource path must be a confined relative path")
    path = (root / relative).resolve(strict=True)
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError("resource path escapes bundle or is not a file")
    raw = path.read_bytes()
    if len(raw) != expected_size or hashlib.sha256(raw).hexdigest() != expected_hash:
        raise ValueError(f"resource size or SHA-256 mismatch: {relative}")
    return raw


def validate_input_bundle(bundle):
    """Check the manifest independently of local file availability."""
    _validate_schema(bundle, "model-input-bundle-v1.schema.json", "ModelInputBundle")
    _verify_content_hash(bundle, "ModelInputBundle")
    if type(bundle["schema_version"]) is not int:
        raise ValueError("schema_version must be an integer")
    host, phage = bundle["host_state"], bundle["phage_profile"]
    _verify_content_hash(host, "HostState")
    if _canonical_sha256(phage) != bundle["phage_profile_sha256"]:
        raise ValueError("PhageProfile hash mismatch")
    created = _time(bundle["created_at"])
    if any(_time(t) > created for t in (host["as_of"], host["created_at"], phage["created_at"])):
        raise ValueError("bundle cannot predate its input snapshots")
    context = bundle["assay_context"]
    if context["status"] == "unknown" and (context["conditions"] or context["evidence"]):
        raise ValueError("unknown assay_context must have no declared conditions/evidence")
    if context["status"] == "declared" and (not context["conditions"] or not context["evidence"]):
        raise ValueError("declared assay_context requires conditions and evidence")
    subjects, paths = set(), set()
    for resource in bundle["resources"]:
        subject = resource["subject"]
        if subject in subjects or resource["path"] in paths:
            raise ValueError("duplicate sequence subject or resource path")
        subjects.add(subject)
        paths.add(resource["path"])
        if subject == "host":
            if resource["sequence_sha256"] != host["host_genome_sha256"]:
                raise ValueError("host sequence identity mismatch")
        else:
            sequence = phage["sequence"] or {}
            snapshot = phage["provenance"].get("snapshot_sha256")
            if sequence.get("sha256"):
                if resource["sequence_sha256"] != sequence["sha256"]:
                    raise ValueError("phage sequence identity mismatch")
            elif not (sequence.get("accession") and snapshot):
                raise ValueError("phage sequence requires exact digest or accession/snapshot")
            if sequence.get("accession") and resource["accession"] != sequence["accession"]:
                raise ValueError("phage accession mismatch")
            if snapshot and resource["source_snapshot_sha256"] != snapshot:
                raise ValueError("phage source snapshot mismatch")
    for subject in ("host", "phage"):
        expected = "available" if subject in subjects else "not_provided"
        if bundle["sequence_status"][subject] != expected:
            raise ValueError("sequence_status does not match packaged resources")
    feature_ids = set()
    for entry in bundle["feature_artifacts"]:
        artifact = entry["manifest"]
        _verify_content_hash(artifact, "FeatureArtifact")
        analysis = artifact["analysis"]
        if _canonical_sha256(analysis["parameters"]) != analysis["parameters_sha256"]:
            raise ValueError("FeatureArtifact parameter hash mismatch")
        started, finished = _time(analysis["started_at"]), _time(analysis["finished_at"])
        if not started <= finished <= _time(artifact["created_at"]) <= created:
            raise ValueError("FeatureArtifact chronology mismatch")
        if artifact["artifact_id"] in feature_ids:
            raise ValueError("duplicate FeatureArtifact ID")
        feature_ids.add(artifact["artifact_id"])
        subject = artifact["subject"]
        if subject["kind"] == "host_genome":
            matches = subject["id"] in {host["record_id"], host["host_profile_id"]}
            expected = host["host_genome_sha256"]
        else:
            matches = subject["id"] == phage["profile_id"]
            expected = (phage["sequence"] or {}).get("sha256")
            if expected is None:
                expected = next(
                    (r["sequence_sha256"] for r in bundle["resources"] if r["subject"] == "phage"),
                    None,
                )
        if not matches or subject["sequence_sha256"] != expected:
            raise ValueError("FeatureArtifact subject does not match exact pair")
        for output in artifact["outputs"]:
            relative = f"{entry['directory']}/{output['path']}"
            if relative in paths:
                raise ValueError("duplicate feature resource path")
            paths.add(relative)
    return bundle


def load_input_bundle(path, *, required_sequences=()):
    """Verify all payloads before model-specific preprocessing reads an export."""
    path = Path(path)
    bundle = validate_input_bundle(read_json_object(path))
    if set(required_sequences) - {"host", "phage"}:
        raise ValueError("required sequences must be host or phage")
    for subject in required_sequences:
        if bundle["sequence_status"][subject] != "available":
            raise ValueError(f"required {subject} sequence not provided")
    for resource in bundle["resources"]:
        _file(path.parent, resource["path"], resource["sha256"], resource["size_bytes"])
    for entry in bundle["feature_artifacts"]:
        for output in entry["manifest"]["outputs"]:
            relative = f"{entry['directory']}/{output['path']}"
            _file(path.parent, relative, output["sha256"], output["size_bytes"])
    return bundle


def read_sequence_bytes(bundle_path, subject):
    """Return the verified sequence file bytes for one subject."""
    bundle = load_input_bundle(bundle_path, required_sequences=(subject,))
    resource = next(r for r in bundle["resources"] if r["subject"] == subject)
    return _file(Path(bundle_path).parent, resource["path"], resource["sha256"], resource["size_bytes"])


def validate_pair_prediction(result, bundle):
    validate_input_bundle(bundle)
    _validate_schema(result, "pair-prediction-v1.schema.json", "PairPrediction")
    _verify_content_hash(result, "PairPrediction")
    if type(result["schema_version"]) is not int:
        raise ValueError("schema_version must be an integer")
    bound = (result["bundle_id"], result["bundle_sha256"])
    if bound != (bundle["bundle_id"], bundle["content_sha256"]):
        raise ValueError("prediction does not bind the exact bundle")
    if _time(result["created_at"]) < _time(bundle["created_at"]):
        raise ValueError("prediction cannot predate bundle")
    estimate, producer = result["estimate"], result["producer"]
    if estimate["status"] == "not_estimated":
        if estimate["value"] is not None or estimate["uncertainty"] is not None:
            raise ValueError("not_estimated must have null values and uncertainty")
    else:
        if (
            estimate["value"] is None
            or not producer["artifact_sha256"]
            or not producer["preprocessing_sha256"]
        ):
            raise ValueError("estimated output requires value, model and preprocessing digests")
        interval = estimate["uncertainty"]
        if interval and not interval["lower"] <= estimate["value"] <= interval["upper"]:
            raise ValueError("uncertainty interval must contain value")
        if estimate["status"] == "probability" and (
                not interval or not producer["calibration_sha256"]
                or bundle["endpoint"]["assay_definition_status"] != "verified"):
            raise ValueError(
                "probability requires verified assay definition, calibration digest and uncertainty"
            )
    return result


def build_pair_prediction(
    bundle, *, prediction_id, estimate=None, producer=None, created_at=None, limitations=None
):
    """Record explicit model output, or abstain when no inference is supplied."""
    result = {
        "schema_version": 1, "kind": "pair_prediction", "prediction_id": prediction_id,
        "bundle_id": bundle["bundle_id"], "bundle_sha256": bundle["content_sha256"],
        "created_at": created_at or datetime.now(timezone.utc).isoformat(),
        "estimate": estimate if estimate is not None else {
            "status": "not_estimated", "value": None, "uncertainty": None,
            "reason": "inference_not_available",
        },
        "producer": producer if producer is not None else {
            "name": "unresolved-no-inference", "version": "1", "artifact_sha256": None,
            "preprocessing_sha256": None, "calibration_sha256": None,
        },
        "limitations": limitations if limitations is not None else [
            "Source-endpoint research output; no mechanistic target mapping or cocktail effect.",
        ],
    }
    result["content_sha256"] = canonical_content_sha256(result)
    return validate_pair_prediction(result, bundle)
