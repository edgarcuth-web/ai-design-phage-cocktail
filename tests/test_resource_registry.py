"""Tests for exact local resource resolution; fixture records are software-only."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from phage_cocktail.model_interface import validate_model_request
from phage_cocktail.resource_registry import (
    load_resource_registry,
    resolve_request_records,
)

ROOT = Path(__file__).resolve().parents[1]
REQUEST_FIXTURE = ROOT / "tests" / "fixtures" / "model_request_v1.json"
TIMESTAMP = "2026-10-17T11:00:00+00:00"


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def content_hash(document):
    content = dict(document)
    content.pop("content_sha256", None)
    canonical = json.dumps(
        content, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return path


def reseal_request(request):
    request["content_sha256"] = content_hash(request)
    return validate_model_request(request)


def make_local_registry(root, *, host_status="available", phage_status="available"):
    """Create matching records with unmistakably software-only test provenance."""
    root = Path(root)
    request = load_json(REQUEST_FIXTURE)

    host_ref = request["host_state"]
    host_state = {
        "schema_version": 1,
        "record_id": host_ref["record_id"],
        "host_state_id": host_ref["record_id"].rsplit(":v", 1)[0],
        "revision": 1,
        "host_profile_id": host_ref["host_profile_id"],
        "host_genome_sha256": host_ref["host_genome_sha256"],
        "as_of": host_ref["as_of"],
        "components": [],
        "evidence": [],
        "limitations": ["Software-only test record; not biological evidence."],
        "previous_state": None,
        "created_at": host_ref["as_of"],
    }
    host_state["content_sha256"] = content_hash(host_state)

    phage_profile = {
        "schema_version": 1,
        "profile_id": request["phage_profile_id"],
        "phage_id": request["phage_profile_id"].rsplit(":v", 1)[0],
        "revision": 1,
        "created_at": TIMESTAMP,
        "provenance": {
            "kind": "retrieved",
            "source_db": "software-test-fixture-only",
            "accession": "fixture-not-a-biological-accession",
            "retrieved_at": TIMESTAMP,
            "source_uri": "https://example.invalid/software-test-fixture",
        },
        "sequence": None,
        "source_annotations": {"fixture_only": True},
    }

    request["host_state"]["content_sha256"] = host_state["content_sha256"]
    request["phage_profile_sha256"] = content_hash(phage_profile)
    reseal_request(request)

    host_path = write_json(root / "records" / "host-state.json", host_state)
    phage_path = write_json(root / "records" / "phage-profile.json", phage_profile)
    registry_path = root / "registry.v1.json"
    manifest = {
        "schema_version": 1,
        "registry_id": "software-test-registry-v1",
        "entries": [
            {
                "kind": "host_state",
                "resource_id": host_state["record_id"],
                "content_sha256": host_state["content_sha256"],
                "file_sha256": file_hash(host_path),
                "path": "records/host-state.json",
                "status": host_status,
            },
            {
                "kind": "phage_profile",
                "resource_id": phage_profile["profile_id"],
                "content_sha256": content_hash(phage_profile),
                "file_sha256": file_hash(phage_path),
                "path": "records/phage-profile.json",
                "status": phage_status,
            },
        ],
    }
    write_json(registry_path, manifest)
    return request, registry_path, host_path, phage_path


class ResourceRegistryTests(unittest.TestCase):
    def test_resolves_exact_records_and_returns_audit_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, host_path, phage_path = make_local_registry(directory)
            resolved = resolve_request_records(request, registry_path)

            self.assertEqual(resolved.host_state["record_id"], request["host_state"]["record_id"])
            self.assertEqual(resolved.phage_profile["profile_id"], request["phage_profile_id"])
            self.assertEqual(resolved.host_state_path, host_path.resolve())
            self.assertEqual(resolved.phage_profile_path, phage_path.resolve())
            summary = resolved.summary()
            self.assertEqual(summary["request_sha256"], request["content_sha256"])
            self.assertEqual(summary["registry"]["registry_id"], "software-test-registry-v1")
            self.assertFalse(summary["sequence_assets_resolved"])
            self.assertFalse(summary["biological_validity_established"])
            self.assertEqual(
                summary["resources"]["host_state"]["content_sha256"],
                request["host_state"]["content_sha256"],
            )

    def test_request_hash_must_match_registered_host_state(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, _, _ = make_local_registry(directory)
            request["host_state"]["content_sha256"] = "e" * 64
            reseal_request(request)
            with self.assertRaisesRegex(ValueError, "different content SHA-256"):
                resolve_request_records(request, registry_path)

    def test_host_state_identity_must_match_request_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, _, _ = make_local_registry(directory)
            request["host_state"]["host_profile_id"] = "wrong-host-profile:v1"
            reseal_request(request)
            with self.assertRaisesRegex(ValueError, "host_profile_id does not match"):
                resolve_request_records(request, registry_path)

    def test_raw_file_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, _, phage_path = make_local_registry(directory)
            phage_path.write_text(
                phage_path.read_text(encoding="utf-8").replace(
                    "fixture_only", "modified_fixture"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "file SHA-256 mismatch"):
                resolve_request_records(request, registry_path)

    def test_canonical_content_tampering_is_rejected_even_if_file_hash_is_updated(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, _, phage_path = make_local_registry(directory)
            phage = load_json(phage_path)
            phage["source_annotations"]["fixture_only"] = False
            write_json(phage_path, phage)

            manifest = load_json(registry_path)
            manifest["entries"][1]["file_sha256"] = file_hash(phage_path)
            write_json(registry_path, manifest)

            with self.assertRaisesRegex(ValueError, "canonical content SHA-256 mismatch"):
                resolve_request_records(request, registry_path)

    def test_withdrawn_resource_is_not_resolved(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, _, _ = make_local_registry(
                directory, phage_status="withdrawn"
            )
            with self.assertRaisesRegex(ValueError, "not available"):
                resolve_request_records(request, registry_path)

    def test_path_traversal_is_rejected_before_file_access(self):
        with tempfile.TemporaryDirectory() as directory:
            _, registry_path, _, _ = make_local_registry(directory)
            manifest = load_json(registry_path)
            manifest["entries"][0]["path"] = "../outside.json"
            write_json(registry_path, manifest)
            with self.assertRaisesRegex(ValueError, "must not contain"):
                load_resource_registry(registry_path)

    def test_duplicate_resource_identity_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            _, registry_path, _, _ = make_local_registry(directory)
            manifest = load_json(registry_path)
            manifest["entries"].append(copy.deepcopy(manifest["entries"][0]))
            write_json(registry_path, manifest)
            with self.assertRaisesRegex(ValueError, "duplicate registry entry"):
                load_resource_registry(registry_path)

    def test_cli_resolves_records_and_reports_unresolved_sequence_assets(self):
        with tempfile.TemporaryDirectory() as directory:
            request, registry_path, _, _ = make_local_registry(directory)
            request_path = write_json(Path(directory) / "request.json", request)
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "phage_cocktail.model_interface",
                    "resolve-request-records",
                    str(request_path),
                    "--registry",
                    str(registry_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            summary = json.loads(completed.stdout)
            self.assertEqual(summary["request_id"], request["request_id"])
            self.assertFalse(summary["sequence_assets_resolved"])
            self.assertFalse(summary["biological_validity_established"])

    def test_cli_fails_closed_when_registry_is_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            request, _, _, _ = make_local_registry(directory)
            request_path = write_json(Path(directory) / "request.json", request)
            missing_registry = Path(directory) / "missing-registry.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "phage_cocktail.model_interface",
                    "resolve-request-records",
                    str(request_path),
                    "--registry",
                    str(missing_registry),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("resource registry not found", completed.stderr)


if __name__ == "__main__":
    unittest.main()
