"""Consumer conformance fixtures are software examples, not biological evidence."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import unittest
import tempfile
import shutil

from phage_cocktail.artifacts import canonical_content_sha256
from phage_cocktail.input_bundle import load_input_bundle, build_pair_prediction, validate_pair_prediction
from phage_cocktail.model_interface import validate_model_result

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/input_bundle/bundle.json"


class InputBundleTests(unittest.TestCase):
    def test_reads_frozen_values_and_explicit_missingness(self):
        bundle = load_input_bundle(FIXTURE, required_sequences=("host",))
        self.assertEqual(bundle["sequence_status"]["host"], "available")
        self.assertEqual(bundle["endpoint"]["assay_definition_status"], "unverified")
        with self.assertRaisesRegex(ValueError, "required.*phage"):
            load_input_bundle(FIXTURE, required_sequences=("phage",))

    def test_unresolved_result_binds_exact_bundle(self):
        bundle = load_input_bundle(FIXTURE)
        result = build_pair_prediction(bundle, prediction_id="fixture", created_at="2026-10-19T12:00:00Z")
        self.assertIsNone(result["estimate"]["value"])
        self.assertEqual(result["bundle_sha256"], bundle["content_sha256"])
        result["bundle_sha256"] = "a" * 64
        result["content_sha256"] = canonical_content_sha256(result)
        with self.assertRaisesRegex(ValueError, "exact bundle"):
            validate_pair_prediction(result, bundle)

    def test_uncalibrated_score_is_not_a_probability(self):
        bundle = load_input_bundle(FIXTURE)
        producer = {"name": "fixture", "version": "1", "artifact_sha256": "a" * 64,
                    "preprocessing_sha256": "b" * 64, "calibration_sha256": None}
        estimate = {"status": "uncalibrated_score", "value": 0.7, "uncertainty": None, "reason": "software fixture"}
        result = build_pair_prediction(bundle, prediction_id="score", estimate=estimate, producer=producer,
                                       created_at="2026-10-19T12:00:00Z")
        result["estimate"]["status"] = "probability"
        result["content_sha256"] = canonical_content_sha256(result)
        with self.assertRaisesRegex(ValueError, "probability requires"):
            validate_pair_prediction(result, bundle)

    def test_packaged_schemas_match_repository_copies(self):
        for name in ("model-input-bundle-v1.schema.json", "pair-prediction-v1.schema.json"):
            self.assertEqual((ROOT / "schemas" / name).read_bytes(),
                             (ROOT / "src/phage_cocktail/schemas" / name).read_bytes())

    def test_source_probability_requires_each_gate_and_stays_source_defined(self):
        bundle = load_input_bundle(FIXTURE)
        bundle["endpoint"]["assay_definition_status"] = "verified"  # Software fixture only.
        bundle["content_sha256"] = canonical_content_sha256(bundle)
        estimate = {"status": "probability", "value": 0.7, "reason": "Software fixture only.",
                    "uncertainty": {"lower": 0.6, "upper": 0.8, "method": "fixture"}}
        producer = {"name": "fixture", "version": "1", "artifact_sha256": "a" * 64,
                    "preprocessing_sha256": "b" * 64, "calibration_sha256": "c" * 64}
        result = build_pair_prediction(
            bundle, prediction_id="source-probability", estimate=estimate, producer=producer,
            created_at="2026-10-19T12:00:00Z",
        )
        with self.assertRaisesRegex(ValueError, "ModelResult"):
            validate_model_result(result)
        for missing in ("verified_assay", "uncertainty", "calibration"):
            with self.subTest(missing=missing):
                altered_bundle, altered_result = copy.deepcopy(bundle), copy.deepcopy(result)
                if missing == "verified_assay":
                    altered_bundle["endpoint"]["assay_definition_status"] = "unverified"
                    altered_bundle["content_sha256"] = canonical_content_sha256(altered_bundle)
                    altered_result["bundle_sha256"] = altered_bundle["content_sha256"]
                elif missing == "uncertainty":
                    altered_result["estimate"]["uncertainty"] = None
                else:
                    altered_result["producer"]["calibration_sha256"] = None
                altered_result["content_sha256"] = canonical_content_sha256(altered_result)
                with self.assertRaisesRegex(ValueError, "probability requires"):
                    validate_pair_prediction(altered_result, altered_bundle)

    def test_resealed_endpoint_or_context_cannot_reuse_a_source_result(self):
        bundle = load_input_bundle(FIXTURE)
        result = build_pair_prediction(bundle, prediction_id="bound-source",
                                       created_at="2026-10-19T12:00:00Z")
        for field, value in (
            ("endpoint", {**bundle["endpoint"], "endpoint_id": "different-source:v2"}),
            ("assay_context", {"status": "declared", "conditions": {"fixture": True},
                               "evidence": ["software-fixture"]}),
        ):
            with self.subTest(field=field):
                changed = copy.deepcopy(bundle)
                changed[field] = value
                changed["content_sha256"] = canonical_content_sha256(changed)
                with self.assertRaisesRegex(ValueError, "exact bundle"):
                    validate_pair_prediction(result, changed)

    def test_cli_checks_required_inputs(self):
        result = subprocess.run([sys.executable, "-m", "phage_cocktail.model_interface",
                                 "validate-bundle", str(FIXTURE), "--require-sequence", "host"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_changed_payload_and_resealed_identity_are_rejected(self):
        for mutate, message in (
            (lambda b: b["resources"][0].update(sequence_sha256="e" * 64), "identity"),
            (lambda b: b["resources"][0].update(path="../outside.fasta"), "relative"),
            (lambda b: b["resources"].append(copy.deepcopy(b["resources"][0])), "duplicate"),
        ):
            with self.subTest(message=message), tempfile.TemporaryDirectory() as directory:
                root = Path(directory) / "input"
                shutil.copytree(FIXTURE.parent, root)
                path = root / "bundle.json"
                bundle = json.loads(path.read_text())
                mutate(bundle)
                bundle["content_sha256"] = canonical_content_sha256(bundle)
                path.write_text(json.dumps(bundle))
                with self.assertRaisesRegex(ValueError, message):
                    load_input_bundle(path)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "input"
            shutil.copytree(FIXTURE.parent, root)
            bundle = load_input_bundle(root / "bundle.json")
            entry = bundle["feature_artifacts"][0]
            (root / entry["directory"] / entry["manifest"]["outputs"][0]["path"]).write_text("changed")
            with self.assertRaisesRegex(ValueError, "SHA-256|size"):
                load_input_bundle(root / "bundle.json")
