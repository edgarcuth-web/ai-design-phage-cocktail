"""Validate the shared v1 model exchange JSON Schemas and fixtures."""

import copy
import hashlib
import json
from pathlib import Path
import unittest

from jsonschema import Draft202012Validator, FormatChecker

ROOT = Path(__file__).resolve().parents[1]
SCHEMA_DIR = ROOT / "schemas"
FIXTURE_DIR = ROOT / "tests" / "fixtures"


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def content_hash(document):
    content = dict(document)
    content.pop("content_sha256")
    canonical = json.dumps(
        content, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


class ModelWireSchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.request_schema = load_json(SCHEMA_DIR / "model-request-v1.schema.json")
        cls.result_schema = load_json(SCHEMA_DIR / "model-result-v1.schema.json")
        cls.request = load_json(FIXTURE_DIR / "model_request_v1.json")
        cls.result = load_json(FIXTURE_DIR / "model_result_v1.json")
        cls.request_validator = Draft202012Validator(
            cls.request_schema, format_checker=FormatChecker()
        )
        cls.result_validator = Draft202012Validator(
            cls.result_schema, format_checker=FormatChecker()
        )

    def test_schemas_are_valid_draft_2020_12(self):
        Draft202012Validator.check_schema(self.request_schema)
        Draft202012Validator.check_schema(self.result_schema)

    def test_shared_fixtures_validate(self):
        self.request_validator.validate(self.request)
        self.result_validator.validate(self.result)

    def test_fixtures_preserve_content_hash_and_request_lineage(self):
        self.assertEqual(self.request["content_sha256"], content_hash(self.request))
        self.assertEqual(self.result["content_sha256"], content_hash(self.result))
        self.assertEqual(self.request["record_id"], self.request["request_id"])
        self.assertEqual(self.result["record_id"], self.result["result_id"])
        self.assertEqual(self.result["request_id"], self.request["request_id"])
        self.assertEqual(
            self.result["request_sha256"], self.request["content_sha256"]
        )
        self.assertEqual(self.result["host_state"], self.request["host_state"])
        self.assertEqual(
            self.result["phage_profile_sha256"],
            self.request["phage_profile_sha256"],
        )

    def test_request_rejects_duplicate_or_unknown_targets(self):
        duplicate = copy.deepcopy(self.request)
        duplicate["requested_targets"].append("productive_infection")
        self.assertFalse(self.request_validator.is_valid(duplicate))

        unknown = copy.deepcopy(self.request)
        unknown["requested_targets"] = ["invented_target"]
        self.assertFalse(self.request_validator.is_valid(unknown))

    def test_unknown_properties_are_rejected(self):
        request = copy.deepcopy(self.request)
        request["unexpected"] = True
        self.assertFalse(self.request_validator.is_valid(request))

        result = copy.deepcopy(self.result)
        result["unexpected"] = True
        self.assertFalse(self.result_validator.is_valid(result))

    def test_not_estimated_target_cannot_carry_placeholder_probability(self):
        result = copy.deepcopy(self.result)
        result["probabilities"]["adsorption"]["probability"] = 0.0
        self.assertFalse(self.result_validator.is_valid(result))

    def test_estimated_target_requires_uncertainty_and_model_artifact(self):
        no_interval = copy.deepcopy(self.result)
        no_interval["probabilities"]["productive_infection"]["uncertainty_lower"] = None
        self.assertFalse(self.result_validator.is_valid(no_interval))

        no_artifact = copy.deepcopy(self.result)
        no_artifact["producer"]["artifact_sha256"] = None
        self.assertFalse(self.result_validator.is_valid(no_artifact))

    def test_uncertainty_bounds_must_be_probabilities(self):
        result = copy.deepcopy(self.result)
        result["probabilities"]["productive_infection"]["uncertainty_upper"] = 1.1
        self.assertFalse(self.result_validator.is_valid(result))


if __name__ == "__main__":
    unittest.main()
