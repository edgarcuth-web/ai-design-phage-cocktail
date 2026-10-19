"""Contract tests for the real PADLOC adapter using a controlled CLI stub.

The stub exercises command construction and output handling only. Actual PADLOC
execution is covered by the dedicated PADLOC integration workflow.
"""

import csv
from contextlib import redirect_stdout
import copy
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from phage_cocktail.artifacts import canonical_content_sha256, sha256_file
from phage_cocktail.feature_artifact import load_feature_artifact
from phage_cocktail.model_interface import main
from phage_cocktail.padloc_adapter import (
    PADLOC_DB_HASH_METHOD,
    PADLOC_DB_RELEASE_URI,
    PADLOC_DB_VERSION,
    PADLOC_VERSION,
    load_padloc_summary,
    run_padloc_host_defense,
)
from phage_cocktail.sequence_registry import fingerprint_fasta

TIMESTAMP = "2026-10-18T12:00:00+00:00"
CSV_FIELDS = [
    "system.number",
    "seqid",
    "system",
    "target.name",
    "hmm.accession",
    "hmm.name",
]


def write_json(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return path


def reseal(document):
    document["content_sha256"] = canonical_content_sha256(document)
    return document


def make_database(root, *, release="v2.0.0"):
    root = Path(root)
    (root / "hmm").mkdir(parents=True, exist_ok=True)
    (root / "cm").mkdir(parents=True, exist_ok=True)
    (root / "sys").mkdir(parents=True, exist_ok=True)
    (root / "RELEASE.md").write_text(
        f"# PADLOC-DB release {release}\n",
        encoding="utf-8",
    )
    (root / "hmm" / "padlocdb.hmm").write_text("HMM fixture bytes\n", encoding="utf-8")
    (root / "cm" / "padlocdb.cm").write_text("CM fixture bytes\n", encoding="utf-8")
    (root / "hmm_meta.txt").write_text("hmm metadata fixture\n", encoding="utf-8")
    (root / "sys_meta.txt").write_text("system metadata fixture\n", encoding="utf-8")
    (root / "sys" / "example.yaml").write_text("name: fixture\n", encoding="utf-8")
    return root


def make_fake_executable(root):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    executable = root / "padloc"
    executable.write_text(
        """#!/usr/bin/env python3
import csv
import os
from pathlib import Path
import sys
import time

args = sys.argv[1:]
if "--version" in args:
    print("padloc " + os.environ.get("PADLOC_TEST_VERSION", "v2.0.0"))
    raise SystemExit(0)
if "--db-version" in args:
    print("padloc-db " + os.environ.get("PADLOC_TEST_DB_VERSION", "v2.0.0"))
    raise SystemExit(0)
if "--check-deps" in args:
    if os.environ.get("PADLOC_TEST_MISSING_DEP") == "1":
        print("WARNING: Prodigal not installed")
    else:
        print("R installed (4.3.1)")
        print("HMMER installed (3.3.2)")
        print("Prodigal installed (2.6.3)")
        print("R package 'tidyverse' installed (2.0.0)")
        print("R package 'yaml' installed (2.3.8)")
        print("R package 'getopt' installed (1.20.4)")
    raise SystemExit(0)
if "--fna" in args:
    if os.environ.get("PADLOC_TEST_MODE") == "timeout":
        time.sleep(5)
    fasta = Path(args[args.index("--fna") + 1])
    outdir = Path(args[args.index("--outdir") + 1])
    outdir.mkdir(parents=True, exist_ok=True)
    assert fasta.is_file()
    (outdir / "host.domtblout").write_text("# HMMER domain table fixture\\n")
    (outdir / "host_prodigal.faa").write_text(">gene_1\\nMKK\\n")
    (outdir / "host_prodigal.gff").write_text("##gff-version 3\\n")
    mode = os.environ.get("PADLOC_TEST_MODE", "hits")
    if mode in {"hits", "empty_csv"}:
        rows = [
            ["1", "contig1", "Restriction-modification", "gene_1", "HMM0001", "fixture_hmm_1"],
            ["1", "contig1", "Restriction-modification", "gene_2", "HMM0002", "fixture_hmm_2"],
            ["2", "contig1", "Abortive infection", "gene_3", "HMM0003", "fixture_hmm_3"],
        ] if mode == "hits" else []
        with (outdir / "host_padloc.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["system.number", "seqid", "system", "target.name", "hmm.accession", "hmm.name"])
            writer.writerows(rows)
        if rows:
            (outdir / "host_padloc.gff").write_text("##gff-version 3\\n")
            print("Complete: Writing output to host_padloc.csv")
        elif mode == "empty_csv":
            print("PADLOC completed")
    elif mode == "no_hits":
        print("Complete: Nothing found for host")
    elif mode == "silent":
        print("PADLOC completed without an output result")
    raise SystemExit(0)
print("unexpected test invocation", args, file=sys.stderr)
raise SystemExit(3)
""",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    return executable


def make_host_fixture(root):
    root = Path(root)
    inputs = root / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    fasta_path = inputs / "host.fna"
    fasta_path.write_text(">contig1\nACGTACGTACGT\n", encoding="utf-8")
    asset = fingerprint_fasta(fasta_path, kind="host_genome")
    host_state = {
        "schema_version": 1,
        "record_id": "host-state-padloc-test:v1",
        "host_state_id": "host-state-padloc-test",
        "revision": 1,
        "host_profile_id": "host-profile-padloc-test:v1",
        "host_genome_sha256": asset.sequence_sha256,
        "as_of": TIMESTAMP,
        "components": [],
        "evidence": [],
        "limitations": ["Software-only adapter test fixture; not biological evidence."],
        "previous_state": None,
        "created_at": TIMESTAMP,
    }
    reseal(host_state)
    host_state_path = write_json(inputs / "host-state.json", host_state)
    sequence_registry = {
        "schema_version": 1,
        "registry_id": "padloc-test-sequence-registry-v1",
        "entries": [
            {
                "kind": "host_genome",
                "sequence_sha256": asset.sequence_sha256,
                "file_sha256": asset.file_sha256,
                "path": "host.fna",
                "status": "available",
            }
        ],
    }
    sequence_registry_path = write_json(
        inputs / "sequence-registry.v1.json", sequence_registry
    )
    return {
        "root": root,
        "inputs": inputs,
        "host_state": host_state,
        "host_state_path": host_state_path,
        "sequence_registry_path": sequence_registry_path,
        "fasta_path": fasta_path,
        "asset": asset,
    }


class PadlocAdapterTests(unittest.TestCase):
    def make_environment(self, root):
        root = Path(root)
        fixture = make_host_fixture(root)
        database = make_database(root / "padloc-data")
        executable = make_fake_executable(root / "bin")
        return fixture, database, executable

    def run_adapter(self, fixture, database, executable, output, **kwargs):
        return run_padloc_host_defense(
            host_state_path=fixture["host_state_path"],
            sequence_registry_path=fixture["sequence_registry_path"],
            padloc_data_dir=database,
            output_dir=output,
            padloc_executable=executable,
            cpu=2,
            timeout_seconds=kwargs.pop("timeout_seconds", 60),
            **kwargs,
        )

    def test_actual_cli_contract_builds_verifiable_feature_artifact_for_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            with patch.dict(os.environ, {"PADLOC_TEST_MODE": "hits"}, clear=False):
                result = self.run_adapter(fixture, database, executable, output)

            summary = load_padloc_summary(result["summary_path"])
            artifact = load_feature_artifact(result["manifest_path"])
            self.assertEqual(summary["outcome"], "hits_detected")
            self.assertEqual(summary["systems_detected"], 2)
            self.assertEqual(summary["gene_hit_count"], 3)
            self.assertEqual(
                summary["system_names"],
                ["Abortive infection", "Restriction-modification"],
            )
            self.assertEqual(summary["analyzer"]["version"], PADLOC_VERSION)
            self.assertEqual(summary["analyzer"]["database"]["version"], PADLOC_DB_VERSION)
            self.assertEqual(
                summary["analyzer"]["database"]["source_uri"], PADLOC_DB_RELEASE_URI
            )
            self.assertEqual(
                summary["analyzer"]["database"]["hash_method"], PADLOC_DB_HASH_METHOD
            )
            self.assertEqual(
                artifact["subject"]["source_record_sha256"],
                fixture["host_state"]["content_sha256"],
            )
            self.assertEqual(
                artifact["subject"]["source_record_file_sha256"],
                sha256_file(fixture["host_state_path"]),
            )
            self.assertEqual(artifact["feature_set"], "host_defense_systems")
            self.assertEqual(artifact["status"], "complete")
            self.assertEqual(len(artifact["outputs"]), 11)
            self.assertTrue((output / "outputs" / "padloc-native" / "host.domtblout").is_file())
            self.assertTrue((output / "outputs" / "padloc-native" / "host_padloc.csv").is_file())

    def test_successful_zero_hit_run_is_explicit_not_fabricated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            with patch.dict(os.environ, {"PADLOC_TEST_MODE": "no_hits"}, clear=False):
                result = self.run_adapter(fixture, database, executable, output)
            summary = load_padloc_summary(result["summary_path"])
            load_feature_artifact(result["manifest_path"])
            self.assertEqual(summary["outcome"], "no_hits_detected")
            self.assertEqual(summary["systems_detected"], 0)
            self.assertEqual(summary["gene_hit_count"], 0)
            self.assertEqual(summary["records"], [])
            self.assertIsNone(summary["native_output_path"])
            self.assertFalse((output / "outputs" / "padloc-native" / "host_padloc.csv").exists())

    def test_cli_runs_adapter_and_prints_auditable_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            stdout = io.StringIO()
            with patch.dict(os.environ, {"PADLOC_TEST_MODE": "hits"}, clear=False):
                with redirect_stdout(stdout):
                    code = main(
                        [
                            "annotate-host-defense-padloc",
                            "--host-state", str(fixture["host_state_path"]),
                            "--sequence-registry", str(fixture["sequence_registry_path"]),
                            "--padloc-data", str(database),
                            "--outdir", str(output),
                            "--padloc-executable", str(executable),
                            "--cpu", "2",
                        ]
                    )
            self.assertEqual(code, 0)
            response = json.loads(stdout.getvalue())
            self.assertEqual(response["outcome"], "hits_detected")
            self.assertTrue(Path(response["manifest_path"]).is_file())

    def test_padloc_version_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            with patch.dict(
                os.environ,
                {"PADLOC_TEST_VERSION": "v9.9.9"},
                clear=False,
            ):
                with self.assertRaisesRegex(ValueError, "PADLOC version mismatch"):
                    self.run_adapter(fixture, database, executable, root / "run-output")

    def test_database_version_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            with patch.dict(
                os.environ,
                {"PADLOC_TEST_DB_VERSION": "v1.3.0"},
                clear=False,
            ):
                with self.assertRaisesRegex(ValueError, "PADLOC-DB version mismatch"):
                    self.run_adapter(fixture, database, executable, root / "run-output")

    def test_database_tree_must_match_pinned_release_marker(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            (database / "RELEASE.md").write_text(
                "# PADLOC-DB release v1.3.0\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "release marker"):
                self.run_adapter(fixture, database, executable, root / "run-output")

    def test_missing_runtime_dependency_fails_and_preserves_diagnostic_log(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            with patch.dict(
                os.environ,
                {"PADLOC_TEST_MISSING_DEP": "1"},
                clear=False,
            ):
                with self.assertRaisesRegex(ValueError, "dependency check reported"):
                    self.run_adapter(fixture, database, executable, output)
            self.assertTrue((output / "logs" / "padloc-check-deps.txt").is_file())
            self.assertFalse((output / "feature-artifact.v1.json").exists())

    def test_host_state_hash_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            changed = copy.deepcopy(fixture["host_state"])
            changed["limitations"] = ["tampered"]
            write_json(fixture["host_state_path"], changed)
            with self.assertRaisesRegex(ValueError, "HostState content_sha256 mismatch"):
                self.run_adapter(fixture, database, executable, root / "run-output")

    def test_host_genome_must_match_host_state_fingerprint(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            changed = copy.deepcopy(fixture["host_state"])
            changed["host_genome_sha256"] = "a" * 64
            reseal(changed)
            write_json(fixture["host_state_path"], changed)
            with self.assertRaisesRegex(ValueError, "no registered host_genome sequence"):
                self.run_adapter(fixture, database, executable, root / "run-output")

    def test_missing_native_csv_without_explicit_no_hit_message_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            with patch.dict(os.environ, {"PADLOC_TEST_MODE": "silent"}, clear=False):
                with self.assertRaisesRegex(ValueError, "neither its native CSV nor an explicit no-hit"):
                    self.run_adapter(fixture, database, executable, output)
            self.assertTrue((output / "logs" / "padloc.stdout.txt").is_file())
            self.assertFalse((output / "feature-artifact.v1.json").exists())

    def test_zero_row_csv_without_no_hit_message_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            with patch.dict(os.environ, {"PADLOC_TEST_MODE": "empty_csv"}, clear=False):
                with self.assertRaisesRegex(ValueError, "zero-row CSV"):
                    self.run_adapter(fixture, database, executable, output)

    def test_nonempty_output_directory_is_never_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            output = root / "run-output"
            output.mkdir()
            (output / "old.txt").write_text("old run", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "must be empty"):
                self.run_adapter(fixture, database, executable, output)

    def test_invalid_cpu_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture, database, executable = self.make_environment(root)
            with self.assertRaisesRegex(ValueError, "cpu must be a positive integer"):
                run_padloc_host_defense(
                    host_state_path=fixture["host_state_path"],
                    sequence_registry_path=fixture["sequence_registry_path"],
                    padloc_data_dir=database,
                    output_dir=root / "run-output",
                    padloc_executable=executable,
                    cpu=0,
                )


if __name__ == "__main__":
    unittest.main()
