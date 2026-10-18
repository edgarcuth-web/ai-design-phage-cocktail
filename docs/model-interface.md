# External model exchange v1

The self-driving lab owns the scientific input records and result admission. This
repository owns data, feature extraction, training, evaluation, calibration,
model artifacts and inference. The exchange is versioned JSON; the repositories
do not import each other's private Python modules.

## Machine-readable contract

- [ModelRequest v1 JSON Schema](../schemas/model-request-v1.schema.json)
- [ModelResult v1 JSON Schema](../schemas/model-result-v1.schema.json)
- [Request conformance fixture](../tests/fixtures/model_request_v1.json)
- [Result conformance fixture](../tests/fixtures/model_result_v1.json)

The lab repository mirrors these schema files under its own `schemas/` directory.
The schema `$id` values point to this repository's canonical paths on `main`.
Schema changes must be coordinated across both repositories and keep the shared
fixtures in sync.

## ModelRequest

A request contains:

- a content-addressed HostState reference, including the exact state hash and
  host genome hash;
- the PhageProfile ID and content fingerprint;
- a non-empty, unique list of requested targets;
- request ID, record ID, creation time, schema version and content SHA-256.

The request intentionally carries identifiers and fingerprints, not raw host or
phage sequence data. The model implementation must resolve approved local
resources or fail explicitly; it must not invent features from IDs alone.

## ModelResult

A result must echo the exact request ID and request SHA-256, preserve the same
HostState and PhageProfile lineage, and include exactly three target slots:
`adsorption`, `productive_infection`, and `resistance_emergence`.

Each target is either:

- `not_estimated`, with probability and uncertainty fields set to `null`; or
- `estimated`, with a probability in [0, 1], an uncertainty interval, and a
  non-empty uncertainty-method description.

Any estimated target requires a non-null immutable model artifact SHA-256 in
`producer.artifact_sha256`. Targets that were not requested must remain
`not_estimated`.

## Schema validation versus admission

JSON Schema validates the wire shape and local field constraints. It cannot by
itself prove that content hashes are correct, that `record_id` equals the
corresponding ID, that a result echoes the exact request, or that the requested
targets match the estimated outputs. The lab's semantic validator remains the
authoritative admission gate for those cross-record constraints.

The fixture values are software conformance examples only. They are not
biological predictions and must not be used as training data or scientific
evidence.

## File-based CLI

Install the project, then use the interface CLI:

```bash
pip install -e ".[dev]"

phage-model-interface validate-request tests/fixtures/model_request_v1.json
phage-model-interface validate-result tests/fixtures/model_result_v1.json \
  --request tests/fixtures/model_request_v1.json

phage-model-interface create-unresolved-result \
  tests/fixtures/model_request_v1.json /tmp/model-result-unresolved.json

phage-model-interface validate-result /tmp/model-result-unresolved.json \
  --request tests/fixtures/model_request_v1.json

pytest
```

`validate-request` and `validate-result` check JSON Schema, canonical content
hashes, derived record IDs, timestamp ordering, and (when the request is supplied)
exact request/result lineage and target scope. Output files are written atomically
only after validation.

`create-unresolved-result` is an explicit transport smoke test, not an inference
command. It writes all three targets as `not_estimated`, uses no model artifact
digest, and records why no estimate was made. Missing files, malformed JSON,
unsupported target names, estimated unrequested targets, invalid uncertainty
intervals, and estimated outputs without a model artifact are rejected.

The schemas are included as package data so the installed CLI does not depend on
the repository checkout being present. CI checks the packaged copies against the
canonical schemas in `schemas/`.

No biological feature resolution, model inference, training, or calibration is
implemented by this interface milestone. A request carries identifiers and
fingerprints, not raw sequence data; an actual inference runner must resolve
approved versioned resources and fail explicitly when they are unavailable.
