import hashlib
import json

import pandas as pd
import pytest

from phage_cocktail.registry import load_dataset_manifest, verify_registered_dataset


REVISION = "f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f"


def _manifest_for(path, frame):
    payload = path.read_bytes()
    summary = {
        "interactions": len(frame),
        "bacteria": frame["bacterium"].nunique(),
        "phages": frame["phage"].nunique(),
        "positive": int(frame["infection"].sum()),
        "negative": int(len(frame) - frame["infection"].sum()),
    }
    return {
        "schema_version": 1,
        "dataset_id": "fixture",
        "version": f"upstream-commit:{REVISION}",
        "source": {
            "project": "Fixture",
            "repository_url": "https://example.invalid/source",
            "revision": REVISION,
            "source_path": "data.csv",
            "download_url": "https://example.invalid/data.csv",
            "citation": "Test fixture only.",
        },
        "local_path": "data/raw/fixture.csv",
        "sha256": hashlib.sha256(payload).hexdigest(),
        "license_status": "not_verified",
        "column_mapping": {
            "bacterium": "bacterium",
            "phage": "phage",
            "infection": "infection",
        },
        "label_semantics": {
            "source_column": "infection",
            "positive_value": 1,
            "positive_meaning": "Observed positive label.",
            "negative_value": 0,
            "negative_meaning": "Observed negative label.",
            "operational_definition_status": "not_verified",
            "operational_definition": "Fixture definition not independently verified.",
            "missing_pair_semantics": "unknown_not_negative",
            "model_target_mapping_status": "not_verified",
            "model_target_mapping_note": "No direct mapping is implied.",
        },
        "expected_summary": summary,
    }


def _write_fixture(tmp_path):
    root = tmp_path
    data_path = root / "data" / "raw" / "fixture.csv"
    data_path.parent.mkdir(parents=True)
    frame = pd.DataFrame(
        [
            {"bacterium": "B1", "phage": "P1", "infection": 1},
            {"bacterium": "B1", "phage": "P2", "infection": 0},
            {"bacterium": "B2", "phage": "P1", "infection": 0},
        ]
    )
    frame.to_csv(data_path, index=False)
    manifest_path = root / "data" / "registry" / "fixture.v1.json"
    manifest_path.parent.mkdir(parents=True)
    manifest = _manifest_for(data_path, frame)
    manifest["local_path"] = "data/raw/fixture.csv"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    return manifest_path, data_path, manifest


def test_manifest_records_unverified_label_to_model_target_mapping(tmp_path):
    manifest_path, _, _ = _write_fixture(tmp_path)

    manifest = load_dataset_manifest(manifest_path)

    assert manifest["label_semantics"]["missing_pair_semantics"] == "unknown_not_negative"
    assert manifest["label_semantics"]["model_target_mapping_status"] == "not_verified"


def test_verify_registered_dataset_checks_hash_and_summary(tmp_path):
    manifest_path, data_path, manifest = _write_fixture(tmp_path)

    verified = verify_registered_dataset(manifest_path)

    assert verified.path == data_path
    assert verified.sha256 == manifest["sha256"]
    assert verified.summary["interactions"] == 3
    assert verified.summary["positive"] == 1


def test_verify_registered_dataset_rejects_hash_mismatch(tmp_path):
    manifest_path, data_path, manifest = _write_fixture(tmp_path)
    manifest["sha256"] = "0" * 64
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_registered_dataset(manifest_path)


def test_verify_registered_dataset_rejects_summary_mismatch(tmp_path):
    manifest_path, _, manifest = _write_fixture(tmp_path)
    manifest["expected_summary"]["interactions"] = 99
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="summary mismatch"):
        verify_registered_dataset(manifest_path)


def test_manifest_rejects_implicit_negative_missing_pair_semantics(tmp_path):
    manifest_path, _, manifest = _write_fixture(tmp_path)
    manifest["label_semantics"]["missing_pair_semantics"] = "negative"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="unknown_not_negative"):
        load_dataset_manifest(manifest_path)
