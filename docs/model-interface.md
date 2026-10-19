# Model exchange boundary

This repository owns model-ready datasets, model-specific preprocessing and
representations, training, evaluation, calibration, model artifacts and inference.
It does not own biological characterization, source-record resolution, raw FASTA
ingestion, PADLOC or other scientific annotation tools. Those responsibilities
belong to `self-driving-lab-for-phage-therapy`.

The current ModelRequest v1 / ModelResult v1 interface is transitional. The target
boundary is an immutable, content-addressed `ModelInputBundle v1` produced by SDL.
That bundle should carry admitted scientific annotations, explicit screening and
missingness states, field-level provenance, and exact source fingerprints. This
repository will then perform only model-specific transformation and inference
against that frozen input.

## Wire schemas

- [ModelRequest v1 JSON Schema](../schemas/model-request-v1.schema.json)
- [ModelResult v1 JSON Schema](../schemas/model-result-v1.schema.json)
- [ModelRequest fixture](../tests/fixtures/model_request_v1.json)
- [ModelResult fixture](../tests/fixtures/model_result_v1.json)

The schemas are mirrored in the SDL repository so both sides can validate the
same wire format. Changes must be coordinated. Schema validation checks the wire
shape; the Python interface additionally checks content hashes, record-ID
equality, exact request/result lineage, timestamp ordering and requested-target
admission.

## Current v1 semantics

A ModelRequest binds one exact HostState content hash, one exact PhageProfile
content fingerprint and the prediction targets being requested. A ModelResult
echoes the exact request ID and SHA-256, the same host/phage fingerprints, all
three target slots, producer identity, immutable model artifact SHA-256,
limitations, explanation and prediction time. Targets that were not requested
must remain `not_estimated`.

This version does not carry the complete scientific characterization bundle yet.
Do not infer that a missing annotation is biologically negative. In particular,
PADLOC `no_hits_detected` means no calls were detected by that pinned analysis
run; it is not evidence that the host lacks anti-phage defense systems. Until
ModelInputBundle v1 is implemented, the model must not attempt to resolve local
HostState, PhageProfile or sequence-registry paths.

## CLI

Validate a request:

```bash
phage-model-interface validate-request tests/fixtures/model_request_v1.json
```

Validate a result against its exact request:

```bash
phage-model-interface validate-result tests/fixtures/model_result_v1.json \
  --request tests/fixtures/model_request_v1.json
```

Create an explicit unresolved result without running inference:

```bash
phage-model-interface create-unresolved-result \
  tests/fixtures/model_request_v1.json /tmp/model-result-unresolved.json
```

The unresolved command never fabricates probabilities. All three targets remain
`not_estimated` unless an actual model implementation is integrated and its
outputs pass the same validation contract.

## Model-specific artifacts

Embeddings, tokenizer output, graph representations, normalization parameters,
training tensors, split assignments, model checkpoints, calibration artifacts and
inference outputs belong here because their meaning depends on a specific model or
training pipeline. They must be versioned and tied to the exact input bundle,
preprocessing configuration and model artifact digest.

Scientific annotations such as receptor/surface characterization, PADLOC defense
calls, phage RBP annotations and experimentally observed phenotypes do not belong
here. SDL stores those records with source provenance; this repository consumes a
frozen export rather than reimplementing the biology.
