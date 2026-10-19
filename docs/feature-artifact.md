# FeatureArtifact v1

## Purpose and boundary

FeatureArtifact is the model repository's provenance manifest for the output of
a computational feature-analysis run. It binds a specific host-genome or phage
sequence to a named feature set, exact tool identity, parameters, reference-data
snapshots, output files and explicit limitations.

This is a metadata and integrity contract, not an analyzer. Building or
verifying a FeatureArtifact does not run a tool, interpret its output, certify a
source annotation, or establish biological or experimental validity. No
feature-extraction tool is automatically invoked by this milestone, and no
feature values are bundled as if they had been computed.

The canonical JSON Schema is
[`feature-artifact-v1.schema.json`](../schemas/feature-artifact-v1.schema.json);
the installed package carries a byte-identical copy.

## Manifest fields

Every manifest contains:

- `schema_version`: integer `1`;
- `artifact_id`: caller-owned stable identifier for this analysis output;
- `created_at`: timezone-aware ISO-8601 timestamp;
- `feature_set`: stable feature namespace/name, for example
  `host_defense_systems` or `phage_rbps`;
- `subject`: exact sequence identity and input-file provenance;
- `analysis`: tool, parameters, execution times and reference snapshots;
- `outputs`: one primary output plus optional raw, log or auxiliary files;
- `status`: `complete` or `partial`, describing computational completion only;
- `limitations`: at least one explicit limitation;
- `content_sha256`: canonical digest of the complete manifest excluding this field.

### Subject identity

`subject.kind` is `host_genome` or `phage_sequence`. The subject records the
normalized `sequence_sha256`, raw `input_file_sha256`, base count, contig count
and caller-owned ID. These two digests are not interchangeable: the normalized
sequence digest follows the lab's FASTA convention, while the file digest covers
the exact input bytes including headers, case, wrapping and line endings.

When known, the subject may also retain an accession, source-snapshot digest,
the sequence registry ID plus the exact registry-file digest, and the content
and raw-file digests of the upstream source record (for example, HostState).
The two sequence-registry fields must be present together. The subject metadata
must come from a sequence asset that has already passed the sequence resolver;
this module validates the lineage fields but does not re-parse the input FASTA.

### Analysis provenance

`analysis.tool` requires a tool name and version. It can additionally record
the executable SHA-256 or a container image digest. `analysis.parameters` is a
JSON object and `parameters_sha256` is recomputed from canonical JSON using
sorted keys, compact separators, UTF-8 and no non-standard numeric constants.

`started_at` and `finished_at` must be timezone-aware and ordered; the
manifest's `created_at` cannot predate the run's finish time. `reference_data`
is an array of exact reference snapshots, each with a name, version and SHA-256,
plus an optional source URI. The array may be empty only when no external
reference data was used. A given reference name/version can appear only once
within a manifest.

An executable digest or container digest strengthens reproducibility but is not
mandatory in v1. If neither is recorded, tool name/version and parameters
provide version-level provenance, not a guarantee of bit-for-bit reproduction.
Do not put API keys, passwords, access tokens, or other secrets in parameters or
environment metadata.

### Outputs

Each output contains a role (`primary`, `raw`, `log` or `auxiliary`), a
manifest-relative POSIX path, format, media type, SHA-256 and byte count.
Structured outputs may include `schema_sha256` to pin their own data schema.
Exactly one output must have role `primary`. Output files are stored separately
from the manifest; their content is not copied into the manifest.

The output digest is computed from the exact file bytes. Paths must be relative,
unique, and confined to the manifest directory; absolute paths, traversal
components, and symlinks escaping that directory are rejected. The verifier
checks both the digest and byte count. A zero-byte output can be hash-verified,
but must not be interpreted as an assessed-negative result without
format-specific semantics.

## Integrity semantics

The manifest's `content_sha256` uses the same canonical JSON convention as the
other model-side records, excluding the `content_sha256` field itself. The
parameter digest is separate, so a parameter mutation is reported even when
someone recomputes the outer manifest digest.

SHA-256 checks detect accidental or unauthorized content changes relative to a
trusted manifest; they are not digital signatures and do not establish who ran
the tool or whether the tool is scientifically appropriate. Production use
should retain immutable source snapshots, the exact executable/container where
practical, tool logs, and reviewed parameter choices.

`status = complete` means the recorded computation completed according to its
producer. `status = partial` means the output is incomplete. Neither status
means a feature has been experimentally validated. Every artifact must carry
explicit limitations, and downstream admission must preserve partial/unknown
states rather than converting missing output into a negative finding.

## API and CLI

The builder hashes existing output files; it does not create them or execute
analysis software:

```python
from pathlib import Path

from phage_cocktail.feature_artifact import (
    build_feature_artifact,
    write_feature_artifact,
)

manifest = Path("analysis-run/feature-artifact.v1.json")
artifact = build_feature_artifact(
    manifest_path=manifest,
    artifact_id="run-2026-10-18-defense-001",
    feature_set="host_defense_systems",
    subject={
        "kind": "host_genome",
        "id": "local-isolate-001",
        "sequence_sha256": "<normalized-sequence-sha256>",
        "input_file_sha256": "<raw-fasta-file-sha256>",
        "length_bases": 1,
        "contigs": 1,
    },
    tool={"name": "your-reviewed-tool", "version": "pinned-version"},
    parameters={"mode": "replace-with-real-tool-parameters"},
    reference_data=[],
    outputs=[
        {
            "role": "primary",
            "path": "outputs/features.json",
            "format": "json",
            "media_type": "application/json",
        }
    ],
    started_at="2026-10-18T10:00:00+00:00",
    finished_at="2026-10-18T10:02:00+00:00",
    limitations=["Computational annotation; not experimentally validated."],
)
write_feature_artifact(artifact, manifest)
```

The values above are placeholders showing the API shape, not a real biological
record. In an actual run, replace them with the verified sequence summary, real
tool/version/parameters, pinned reference snapshots, observed execution times,
and output files created by that tool.

Verify an existing manifest and all referenced output files:

```bash
phage-model-interface verify-feature-artifact path/to/feature-artifact.v1.json
```

The CLI prints only the artifact ID, manifest digest, computational status and
output count. It does not print feature content or run inference.

## Next boundary

The first real adapter is [PADLOC host-defense analysis](padloc-adapter.md),
which consumes a registry-verified host genome and pins PADLOC/PADLOC-DB v2.0.0.
The next milestone is a separate, explicit mapping from PADLOC's native output
into HostProfile components. That mapping must preserve source rows, coordinates,
HMM identifiers, scores/coverage, database release and the distinction between
not detected, not searched and unknown; it must not treat an absent call as a
biological negative. Additional tools should be added only when their feature
class, reference release, licensing and output semantics are clear.
