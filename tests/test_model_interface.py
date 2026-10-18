"""End-to-end tests for the external model exchange file interface."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from phage_cocktail.model_interface import (
    PREDICTION_TARGETS,
    build_model_result,
    create_unresolved_result,
    load_model_request,
    load_model_result,
    validate_model_result,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_DIR = ROOT / "tests" / "fixtures"
REQUEST_PATH = FIXTURE_DIR / "model_request_v1.json"
RESULT_PATH = FIXTURE_DIR / "model_result_v1.json"
ARTIFACT_SHA256 = "c" * 64


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def reseal(document):
    """Recompute the fixture hash after an intentional semantic mutation."""
    import hashlib

    content = dict(document)
    content.pop("content_sha256", None)
    canonical = json.dumps(
        content, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    document["content_sha256"] = hashlib.sha256(canonical).hexdigest()
    return document


def estimated(probability=0.7, lower=0.6, upper=0.8):
    return {
        "status": "estimated",
        "probability": probability,
        "uncertainty_lower": lower,
        "uncertainty_upper": upper,
        "uncertainty_method": "software-fixture-only",
    }


def producer(artifact_sha256=ARTIFACT_SHA256):
    return {
        "name": "interface-test-model",
        "version": "test",
        "parameters": {"fixture": True},
        "artifact_sha256": artifact_sha256,
    }


class ModelInterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.request = load_model_request(REQUEST_PATH)
        cls.result = load_model_result(RESULT_PATH, request=cls.request)

    def test_packaged_schemas_match_canonical_repository_schemas(self):
        canonical = ROOT / "schemas"
        packaged = ROOT / "src" / "phage_cocktail" / "schemas"
        for name in (
            "model-request-v1.schema.json",
            "model-result-v1.schema.json",
        ):
            self.assertEqual(
                (packaged / name).read_bytes(),
                (canonical / name).read_bytes(),
                f"packaged schema {name} drifted from canonical schema",
            )

    def test_shared_fixtures_validate_with_exact_lineage(self):
        self.assertEqual(self.result["request_sha256"], self.request["content_sha256"])
        self.assertEqual(self.result["host_state"], self.request["host_state"])
        self.assertEqual(
            self.result["phage_profile_sha256"],
            self.request["phage_profile_sha256"],
        )

    def test_request_hash_tampering_is_rejected(self):
        tampered = copy.deepcopy(self.request)
        tampered["phage_profile_id"] = "changed-phage"
        with self.assertRaisesRegex(ValueError, "content_sha256 mismatch"):
            from phage_cocktail.model_interface import validate_model_request

            validate_model_request(tampered)

    def test_unsupported_target_is_rejected_explicitly(self):
        with self.assertRaisesRegex(ValueError, "unsupported prediction target"):
            build_model_result(
                self.request,
                result_id="unsupported-target-result",
                probabilities={"unknown_target": estimated()},
                producer=producer(),
                mechanistic_explanation="Test only.",
                limitations=["Test fixture only."],
            )

    def test_estimated_unrequested_target_is_rejected(self):
        mutated = copy.deepcopy(self.result)
        mutated["probabilities"]["adsorption"] = estimated(0.4, 0.3, 0.5)
        reseal(mutated)
        with self.assertRaisesRegex(ValueError, "unrequested target"):
            validate_model_result(mutated, request=self.request)

    def test_estimated_result_requires_immutable_model_artifact(self):
        with self.assertRaisesRegex(ValueError, "artifact_sha256"):
            build_model_result(
                self.request,
                result_id="missing-artifact-result",
                probabilities={"productive_infection": estimated()},
                producer=producer(artifact_sha256=None),
                mechanistic_explanation="Test only.",
                limitations=["Test fixture only."],
            )

    def test_uncertainty_interval_must_contain_estimate(self):
        with self.assertRaisesRegex(ValueError, "interval must contain"):
            build_model_result(
                self.request,
                result_id="bad-interval-result",
                probabilities={
                    "productive_infection": estimated(0.7, 0.8, 0.9)
                },
                producer=producer(),
                mechanistic_explanation="Test only.",
                limitations=["Test fixture only."],
                predicted_at="2026-10-17T12:02:00Z",
            )

    def test_builder_fills_unprovided_targets_as_not_estimated(self):
        result = build_model_result(
            self.request,
            result_id="partial-result",
            probabilities={"productive_infection": estimated()},
            producer=producer(),
            mechanistic_explanation="Software fixture only.",
            limitations=["No biological claim."],
            predicted_at="2026-10-17T12:02:00Z",
        )
        self.assertEqual(set(result["probabilities"]), set(PREDICTION_TARGETS))
        self.assertEqual(result["probabilities"]["adsorption"]["status"], "not_estimated")
        self.assertEqual(
            result["probabilities"]["resistance_emergence"]["status"],
            "not_estimated",
        )
        validate_model_result(result, request=self.request)

    def test_unresolved_result_is_explicit_and_content_addressed(self):
        result = create_unresolved_result(
            self.request,
            result_id="explicit-unresolved-result",
            predicted_at="2026-10-17T12:02:00Z",
        )
        self.assertIsNone(result["producer"]["artifact_sha256"])
        for target in PREDICTION_TARGETS:
            self.assertEqual(result["probabilities"][target]["status"], "not_estimated")
        validate_model_result(result, request=self.request)

    def test_cli_validates_request_and_result_files(self):
        commands = [
            [sys.executable, "-m", "phage_cocktail.model_interface", "validate-request", str(REQUEST_PATH)],
            [sys.executable, "-m", "phage_cocktail.model_interface", "validate-result", str(RESULT_PATH), "--request", str(REQUEST_PATH)],
        ]
        for command in commands:
            completed = subprocess.run(
                command, cwd=ROOT, capture_output=True, text=True, check=False
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_cli_writes_an_unresolved_result_without_running_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            request_path = Path(directory) / "request.json"
            result_path = Path(directory) / "result.json"
            request_path.write_text(
                json.dumps(self.request, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "phage_cocktail.model_interface",
                    "create-unresolved-result",
                    str(request_path),
                    str(result_path),
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("no inference was run", completed.stdout)
            written = load_model_result(result_path, request=self.request)
            for target in PREDICTION_TARGETS:
                self.assertEqual(written["probabilities"][target]["status"], "not_estimated")

    def test_cli_reports_missing_request_file_cleanly(self):
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "phage_cocktail.model_interface",
                "validate-request",
                str(ROOT / "does-not-exist.json"),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("cannot read valid JSON", completed.stderr)


if __name__ == "__main__":
    unittest.main()
