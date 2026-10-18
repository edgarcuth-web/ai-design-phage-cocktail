"""Tests for content-addressed FASTA resolution; fixtures are software-only."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from phage_cocktail.model_interface import validate_model_request
from phage_cocktail.sequence_registry import (
    load_sequence_registry,
    resolve_request_sequences,
)

ROOT = Path(__file__).resolve().parents[1]
REQUEST_FIXTURE = ROOT / "tests" / "fixtures" / "model_request_v1.json"
TIMESTAMP = "2026-10-17T11:00:00+00:00"
HOST_FASTA = b">host contig 1\nacgt\n>plasmid\nnry\n"
PHAGE_FASTA = b">phage sequence\nacgtac\n"
HOST_SEQUENCE_SHA256 = hashlib.sha256(b"\0ACGT\0NRY").hexdigest()
PHAGE_SEQUENCE_SHA256 = hashlib.sha256(b"\0ACGTAC").hexdigest()


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


def make_fixture(
    root,
    *,
    phage_has_sha=True,
    phage_length=6,
    phage_status="available",
    host_status="available",
):
    """Create a self-consistent model request and software-only local registries."""
    root = Path(root)
    request = load_json(REQUEST_FIXTURE)
    request["phage_profile_id"] = "phage-profile-fixture:v1"
    request["host_state"]["host_genome_sha256"] = HOST_SEQUENCE_SHA256

    host_state = {
        "schema_version": 1,
        "record_id": request["host_state"]["record_id"],
        "host_state_id": "host-state-fixture",
        "revision": 1,
        "host_profile_id": request["host_state"]["host_profile_id"],
        "host_genome_sha256": HOST_SEQUENCE_SHA256,
        "as_of": request["host_state"]["as_of"],
        "components": [],
        "evidence": [],
        "limitations": ["Software-only test record; not biological evidence."],
        "previous_state": None,
        "created_at": request["host_state"]["as_of"],
    }
    host_state["content_sha256"] = content_hash(host_state)
    request["host_state"]["content_sha256"] = host_state["content_sha256"]

    sequence_ref = {"accession": "PHAGE-TEST.1", "length": phage_length}
    if phage_has_sha:
        sequence_ref["sha256"] = PHAGE_SEQUENCE_SHA256
    phage_profile = {
        "schema_version": 1,
        "profile_id": request["phage_profile_id"],
        "phage_id": "phage-profile-fixture",
        "revision": 1,
        "created_at": TIMESTAMP,
        "provenance": {
            "kind": "retrieved",
            "source_db": "software-test-fixture-only",
            "accession": "PHAGE-TEST.1",
            "retrieved_at": TIMESTAMP,
            "source_uri": "https://example.invalid/software-test-fixture",
        },
        "sequence": sequence_ref,
        "source_annotations": {"fixture_only": True},
    }
    request["phage_profile_sha256"] = content_hash(phage_profile)
    reseal_request(request)

    host_record_path = write_json(root / "records" / "host-state.json", host_state)
    phage_record_path = write_json(
        root / "records" / "phage-profile.json", phage_profile
    )
    resource_registry_path = root / "registry.v1.json"
    write_json(
        resource_registry_path,
        {
            "schema_version": 1,
            "registry_id": "software-test-resource-registry-v1",
            "entries": [
                {
                    "kind": "host_state",
                    "resource_id": host_state["record_id"],
                    "content_sha256": host_state["content_sha256"],
                    "file_sha256": file_hash(host_record_path),
                    "path": "records/host-state.json",
                    "status": "available",
                },
                {
                    "kind": "phage_profile",
                    "resource_id": phage_profile["profile_id"],
                    "content_sha256": content_hash(phage_profile),
                    "file_sha256": file_hash(phage_record_path),
                    "path": "records/phage-profile.json",
                    "status": "available",
                },
            ],
        },
    )

    host_fasta_path = root / "sequences" / "host-genome.fna"
    phage_fasta_path = root / "sequences" / "phage.fna"
    host_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    host_fasta_path.write_bytes(HOST_FASTA)
    phage_fasta_path.write_bytes(PHAGE_FASTA)

    sequence_registry_path = root / "sequence-registry.v1.json"
    sequence_entries = [
        {
            "kind": "host_genome",
            "sequence_sha256": HOST_SEQUENCE_SHA256,
            "file_sha256": file_hash(host_fasta_path),
            "path": "sequences/host-genome.fna",
            "status": host_status,
        },
        {
            "kind": "phage_sequence",
            "sequence_sha256": PHAGE_SEQUENCE_SHA256,
            "file_sha256": file_hash(phage_fasta_path),
            "path": "sequences/phage.fna",
            "status": phage_status,
            "accession": "PHAGE-TEST.1",
        },
    ]
    write_json(
        sequence_registry_path,
        {
            "schema_version": 1,
            "registry_id": "software-test-sequence-registry-v1",
            "entries": sequence_entries,
        },
    )
    return {
        "request": request,
        "resource_registry": resource_registry_path,
        "sequence_registry": sequence_registry_path,
        "host_fasta": host_fasta_path,
        "phage_fasta": phage_fasta_path,
        "host_record": host_record_path,
        "phage_record": phage_record_path,
    }


class SequenceRegistryTests(unittest.TestCase):
    def test_resolves_both_sequences_and_returns_audit_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            resolved = resolve_request_sequences(
                fixture["request"],
                fixture["resource_registry"],
                fixture["sequence_registry"],
            )
            self.assertEqual(
                resolved.host_genome.sequence_sha256, HOST_SEQUENCE_SHA256
            )
            self.assertEqual(
                resolved.phage_sequence.sequence_sha256, PHAGE_SEQUENCE_SHA256
            )
            self.assertEqual(resolved.host_genome.bases, 7)
            self.assertEqual(resolved.host_genome.contigs, 2)
            self.assertEqual(resolved.phage_sequence.bases, 6)
            summary = resolved.summary()
            self.assertTrue(summary["sequence_assets_resolved"])
            self.assertFalse(summary["features_derived"])
            self.assertFalse(summary["biological_validity_established"])
            self.assertEqual(
                summary["sequences"]["host_genome"]["file_sha256"],
                file_hash(fixture["host_fasta"]),
            )
            self.assertEqual(
                summary["sequences"]["phage_sequence"]["accession"],
                "PHAGE-TEST.1",
            )

    def test_headers_case_wrapping_and_line_endings_do_not_change_sequence_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            fixture["host_fasta"].write_bytes(
                b">renamed header\r\nAC\r\nGT\r\n>another name\r\nNRY\r\n"
            )
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][0]["file_sha256"] = file_hash(fixture["host_fasta"])
            write_json(fixture["sequence_registry"], manifest)

            resolved = resolve_request_sequences(
                fixture["request"],
                fixture["resource_registry"],
                fixture["sequence_registry"],
            )
            self.assertEqual(
                resolved.host_genome.sequence_sha256, HOST_SEQUENCE_SHA256
            )
            self.assertNotEqual(
                resolved.host_genome.file_sha256, hashlib.sha256(HOST_FASTA).hexdigest()
            )

    def test_raw_file_tampering_is_rejected_even_if_normalized_sequence_is_unchanged(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            fixture["phage_fasta"].write_bytes(b">changed header\nACGTAC\n")
            with self.assertRaisesRegex(ValueError, "file SHA-256 mismatch"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_sequence_change_is_rejected_even_if_file_hash_is_re_registered(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            fixture["phage_fasta"].write_bytes(b">phage sequence\nACGTAT\n")
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][1]["file_sha256"] = file_hash(fixture["phage_fasta"])
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "normalized sequence SHA-256 mismatch"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_accession_must_match_profile_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][1]["accession"] = "OTHER-PHAGE.1"
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "accession does not match"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_profile_length_must_match_resolved_sequence(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory, phage_length=7)
            with self.assertRaisesRegex(ValueError, "sequence length mismatch"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_accession_only_phage_reference_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory, phage_has_sha=False)
            with self.assertRaisesRegex(ValueError, "must include a normalized SHA-256"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_withdrawn_sequence_is_not_resolved(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory, phage_status="withdrawn")
            with self.assertRaisesRegex(ValueError, "not available"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_withdrawn_host_genome_is_not_resolved(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory, host_status="withdrawn")
            with self.assertRaisesRegex(ValueError, "not available"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_duplicate_sequence_identity_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"].append(copy.deepcopy(manifest["entries"][0]))
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "duplicate sequence registry entry"):
                load_sequence_registry(fixture["sequence_registry"])

    def test_empty_contig_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            fixture["host_fasta"].write_bytes(
                b">first contig\nACGT\n>empty contig\n>third contig\nNRY\n"
            )
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][0]["file_sha256"] = file_hash(fixture["host_fasta"])
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "empty contig"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_path_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][0]["path"] = "../outside.fna"
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "must not contain"):
                load_sequence_registry(fixture["sequence_registry"])

    def test_symlink_cannot_escape_sequence_registry_root(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "registry-root"
            fixture = make_fixture(root)
            outside = Path(directory) / "outside.fna"
            outside.write_bytes(fixture["host_fasta"].read_bytes())
            fixture["host_fasta"].unlink()
            fixture["host_fasta"].symlink_to(outside)

            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][0]["file_sha256"] = file_hash(outside)
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "escapes the sequence registry directory"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_invalid_dna_characters_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            fixture["host_fasta"].write_bytes(b">host\nACGTXYZ\n")
            manifest = load_json(fixture["sequence_registry"])
            manifest["entries"][0]["file_sha256"] = file_hash(fixture["host_fasta"])
            write_json(fixture["sequence_registry"], manifest)
            with self.assertRaisesRegex(ValueError, "invalid DNA characters"):
                resolve_request_sequences(
                    fixture["request"],
                    fixture["resource_registry"],
                    fixture["sequence_registry"],
                )

    def test_cli_resolves_sequences_and_reports_no_features_or_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            request_path = write_json(
                Path(directory) / "request.json", fixture["request"]
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "phage_cocktail.model_interface",
                    "resolve-request-sequences",
                    str(request_path),
                    "--registry",
                    str(fixture["resource_registry"]),
                    "--sequence-registry",
                    str(fixture["sequence_registry"]),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            summary = json.loads(completed.stdout)
            self.assertTrue(summary["sequence_assets_resolved"])
            self.assertFalse(summary["features_derived"])
            self.assertFalse(summary["biological_validity_established"])

    def test_cli_fails_closed_when_sequence_registry_is_missing(self):
        with tempfile.TemporaryDirectory() as directory:
            fixture = make_fixture(directory)
            request_path = write_json(
                Path(directory) / "request.json", fixture["request"]
            )
            missing_registry = Path(directory) / "missing-sequence-registry.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "phage_cocktail.model_interface",
                    "resolve-request-sequences",
                    str(request_path),
                    "--registry",
                    str(fixture["resource_registry"]),
                    "--sequence-registry",
                    str(missing_registry),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn("sequence registry not found", completed.stderr)


if __name__ == "__main__":
    unittest.main()
