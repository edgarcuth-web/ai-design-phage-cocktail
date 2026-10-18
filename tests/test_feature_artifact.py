"""FeatureArtifact v1 contract and output-integrity tests."""

import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from phage_cocktail.artifacts import canonical_content_sha256
from phage_cocktail.feature_artifact import (
    build_feature_artifact,
    load_feature_artifact,
    validate_feature_artifact,
    verify_feature_artifact_outputs,
    write_feature_artifact,
)

ROOT = Path(__file__).resolve().parents[1]
STARTED = "2026-10-18T10:00:00+00:00"
FINISHED = "2026-10-18T10:02:00+00:00"
CREATED = "2026-10-18T10:03:00+00:00"
SEQUENCE_SHA256 = hashlib.sha256(b"\0ACGTACGT").hexdigest()
INPUT_FILE_SHA256 = "b" * 64
TOOL_SHA256 = "c" * 64
REFERENCE_SHA256 = "d" * 64
SOURCE_SNAPSHOT_SHA256 = "e" * 64


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def reseal(document):
    document["content_sha256"] = canonical_content_sha256(document)
    return document


def make_artifact(root, *, status="complete", references=True, outputs=None):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    (root / "outputs").mkdir(exist_ok=True)
    primary = root / "outputs" / "features.json"
    primary.write_text('{"features":[],"fixture_only":true}\n', encoding="utf-8")
    raw = root / "outputs" / "tool-output.gff3"
    raw.write_text("##gff-version 3\n# software test fixture only\n", encoding="utf-8")
    output_specs = outputs or [
        {
            "role": "primary",
            "path": "outputs/features.json",
            "format": "json",
            "media_type": "application/json",
        },
        {
            "role": "raw",
            "path": "outputs/tool-output.gff3",
            "format": "gff3",
            "media_type": "text/gff3",
        },
    ]
    manifest_path = root / "feature-artifact.v1.json"
    artifact = build_feature_artifact(
        manifest_path=manifest_path,
        artifact_id="feature-artifact-test-001",
        feature_set="host_defense_systems",
        subject={
            "kind": "host_genome",
            "id": "host-fixture-001",
            "sequence_sha256": SEQUENCE_SHA256,
            "input_file_sha256": INPUT_FILE_SHA256,
            "length_bases": 8,
            "contigs": 1,
            "sequence_registry_id": "test-sequence-registry",
            "sequence_registry_file_sha256": SOURCE_SNAPSHOT_SHA256,
        },
        tool={
            "name": "example-annotation-tool",
            "version": "1.2.3",
            "executable_sha256": TOOL_SHA256,
        },
        parameters={"mode": "fixture", "threshold": 0.5},
        reference_data=(
            [
                {
                    "name": "example-reference-db",
                    "version": "2026-10",
                    "sha256": REFERENCE_SHA256,
                    "source_uri": "https://example.org/reference-db",
                }
            ]
            if references
            else []
        ),
        outputs=output_specs,
        started_at=STARTED,
        finished_at=FINISHED,
        created_at=CREATED,
        status=status,
        limitations=["Software-only fixture; not a biological result."],
        environment={"runtime": "Python 3.12", "platform": "test"},
    )
    return manifest_path, artifact, primary, raw


class FeatureArtifactTests(unittest.TestCase):
    def test_build_write_load_and_verify_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path, artifact, primary, raw = make_artifact(Path(directory))
            write_feature_artifact(artifact, path)
            loaded = load_feature_artifact(path)

            self.assertEqual(loaded, artifact)
            self.assertEqual(loaded["schema_version"], 1)
            self.assertEqual(loaded["subject"]["sequence_sha256"], SEQUENCE_SHA256)
            self.assertEqual(loaded["analysis"]["parameters_sha256"], hashlib.sha256(
                b'{"mode":"fixture","threshold":0.5}'
            ).hexdigest())
            verified = verify_feature_artifact_outputs(loaded, path)
            self.assertTrue(verified["verified"])
            self.assertEqual(len(verified["outputs"]), 2)
            self.assertEqual(
                verified["outputs"][0]["sha256"],
                hashlib.sha256(primary.read_bytes()).hexdigest(),
            )
            self.assertTrue(raw.is_file())

    def test_reference_database_may_be_explicitly_empty(self):
        with tempfile.TemporaryDirectory() as directory:
            path, artifact, _, _ = make_artifact(Path(directory), references=False)
            self.assertEqual(artifact["analysis"]["reference_data"], [])
            write_feature_artifact(artifact, path)
            self.assertEqual(load_feature_artifact(path)["analysis"]["reference_data"], [])

    def test_partial_status_requires_explicit_limitations(self):
        with tempfile.TemporaryDirectory() as directory:
            _, artifact, _, _ = make_artifact(Path(directory), status="partial")
            self.assertEqual(validate_feature_artifact(artifact)["status"], "partial")

    def test_manifest_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            _, artifact, _, _ = make_artifact(Path(directory))
            tampered = copy.deepcopy(artifact)
            tampered["feature_set"] = "changed_feature_set"
            with self.assertRaisesRegex(ValueError, "content_sha256 mismatch"):
                validate_feature_artifact(tampered)

    def test_parameter_hash_is_checked_independently(self):
        with tempfile.TemporaryDirectory() as directory:
            _, artifact, _, _ = make_artifact(Path(directory))
            tampered = copy.deepcopy(artifact)
            tampered["analysis"]["parameters"]["threshold"] = 0.75
            reseal(tampered)
            with self.assertRaisesRegex(ValueError, "parameters_sha256 mismatch"):
                validate_feature_artifact(tampered)

    def test_output_byte_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, artifact, primary, _ = make_artifact(Path(directory))
            write_feature_artifact(artifact, path)
            primary.write_text('{"features":["changed"]}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "output SHA-256 mismatch"):
                load_feature_artifact(path)

    def test_output_size_is_checked_separately(self):
        with tempfile.TemporaryDirectory() as directory:
            path, artifact, _, _ = make_artifact(Path(directory))
            tampered = copy.deepcopy(artifact)
            tampered["outputs"][0]["size_bytes"] += 1
            reseal(tampered)
            with self.assertRaisesRegex(ValueError, "output size mismatch"):
                verify_feature_artifact_outputs(tampered, path)

    def test_missing_or_multiple_primary_outputs_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            no_primary = [
                {
                    "role": "raw",
                    "path": "outputs/tool-output.gff3",
                    "format": "gff3",
                    "media_type": "text/gff3",
                }
            ]
            with self.assertRaisesRegex(ValueError, "invalid FeatureArtifact"):
                make_artifact(root / "none", outputs=no_primary)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            multiple_primary = [
                {
                    "role": "primary",
                    "path": "outputs/features.json",
                    "format": "json",
                    "media_type": "application/json",
                },
                {
                    "role": "primary",
                    "path": "outputs/tool-output.gff3",
                    "format": "gff3",
                    "media_type": "text/gff3",
                },
            ]
            with self.assertRaisesRegex(ValueError, "invalid FeatureArtifact"):
                make_artifact(root, outputs=multiple_primary)

    def test_duplicate_output_paths_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            outputs = [
                {
                    "role": "primary",
                    "path": "outputs/features.json",
                    "format": "json",
                    "media_type": "application/json",
                },
                {
                    "role": "raw",
                    "path": "outputs/features.json",
                    "format": "json",
                    "media_type": "application/json",
                },
            ]
            with self.assertRaisesRegex(ValueError, "duplicate FeatureArtifact output path"):
                make_artifact(Path(directory), outputs=outputs)

    def test_output_path_traversal_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path, artifact, _, _ = make_artifact(Path(directory))
            tampered = copy.deepcopy(artifact)
            tampered["outputs"][0]["path"] = "../outside.json"
            reseal(tampered)
            with self.assertRaisesRegex(ValueError, "must not be absolute or contain"):
                verify_feature_artifact_outputs(tampered, path)

    def test_symlink_cannot_escape_manifest_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            root = base / "artifact"
            path, artifact, _, _ = make_artifact(root)
            outside = base / "outside.json"
            outside.write_text('{"outside":true}\n', encoding="utf-8")
            link = root / "outputs" / "outside-link.json"
            try:
                link.symlink_to(outside)
            except (OSError, NotImplementedError):
                self.skipTest("symlinks are not available in this environment")
            with self.assertRaisesRegex(ValueError, "escapes the FeatureArtifact directory"):
                build_feature_artifact(
                    manifest_path=path,
                    artifact_id="symlink-test",
                    feature_set="test",
                    subject=artifact["subject"],
                    tool=artifact["analysis"]["tool"],
                    parameters=artifact["analysis"]["parameters"],
                    reference_data=artifact["analysis"]["reference_data"],
                    outputs=[
                        {
                            "role": "primary",
                            "path": "outputs/outside-link.json",
                            "format": "json",
                            "media_type": "application/json",
                        }
                    ],
                    started_at=STARTED,
                    finished_at=FINISHED,
                    created_at=CREATED,
                    limitations=["Software test only."],
                )

    def test_analysis_chronology_is_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            _, artifact, _, _ = make_artifact(Path(directory))
            tampered = copy.deepcopy(artifact)
            tampered["analysis"]["finished_at"] = "2026-10-18T09:59:00+00:00"
            tampered["content_sha256"] = canonical_content_sha256(tampered)
            with self.assertRaisesRegex(ValueError, "cannot predate"):
                validate_feature_artifact(tampered)

    def test_reference_snapshot_duplicates_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            _, artifact, _, _ = make_artifact(Path(directory))
            tampered = copy.deepcopy(artifact)
            tampered["analysis"]["reference_data"].append(
                copy.deepcopy(tampered["analysis"]["reference_data"][0])
            )
            reseal(tampered)
            with self.assertRaisesRegex(ValueError, "duplicate reference snapshot"):
                validate_feature_artifact(tampered)

    def test_cli_verifies_artifact_and_rejects_tampered_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path, artifact, primary, _ = make_artifact(Path(directory))
            write_feature_artifact(artifact, path)
            command = [
                sys.executable,
                "-m",
                "phage_cocktail.model_interface",
                "verify-feature-artifact",
                str(path),
            ]
            good = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(good.returncode, 0, good.stderr)
            self.assertIn("valid FeatureArtifact v1", good.stdout)

            primary.write_text('{"features":["tampered"]}\n', encoding="utf-8")
            bad = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(bad.returncode, 2)
            self.assertIn("output SHA-256 mismatch", bad.stderr)


if __name__ == "__main__":
    unittest.main()
