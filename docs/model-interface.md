# Model exchange boundary

This repository owns model-specific preprocessing, datasets, training, evaluation,
calibration, model artifacts and inference. SDL owns scientific characterization,
source resolution, experimental state and decisions.

Neither side imports the other's private implementation. The exchange is
versioned files with mirrored JSON Schemas and shared software fixtures. The
schemas live in `src/phage_cocktail/schemas/` and are installed with the package;
the shared fixture is `tests/fixtures/input_bundle/`. Embedded FeatureArtifact
definitions mirror SDL's owned v1 schema; no analyzer code moves into this
repository.

See [endpoint semantics](endpoint-semantics.md) for readout distinctions, the
historical versus audited MLC labels, and review requirements before mapping a
source result to a mechanistic target. This adds no new wire fields or mappings.

## Current first-task interface

**ModelInputBundle v1 → PairPrediction v1** supports one existing host–phage pair
and one explicitly source-defined binary endpoint. This is a research interface;
it does not implement inference, inverse design or cocktail-effect prediction.

| Bundle field | Meaning |
| --- | --- |
| `host_state` | Full frozen HostState, components, evidence and missingness, with its own hash |
| `phage_profile` / `phage_profile_sha256` | Full profile, sequence reference, annotations and natural/retrieved/engineered provenance |
| `endpoint` | Named binary source label, dataset version/digest, original label column, assay-definition review status and evidence |
| `assay_context` | `unknown` with empty conditions/evidence, or declared conditions with evidence; never invented defaults |
| `sequence_status` | Explicit `available` or `not_provided` for host and phage |
| `resources` | Explicitly exported local DNA FASTA files, bundle-relative path, byte SHA-256/size, normalized sequence identity and registry provenance |
| `feature_artifacts` | Unchanged scientific manifests plus copied outputs under a bundle-relative directory; partial status, limitations, analyzer/reference versions and output hashes are retained |

All resources live below the export directory. Paths cannot escape through `..`,
absolute paths or symlinks. The consumer verifies every supplied resource and
fails when a model-required sequence is missing. This revision has no network
fetcher or automatic external-service upload. The user must explicitly select
sequences to export; a profile hash alone never implies sequence availability.

`feature_artifacts=[]` means no supplied annotation artifact, not a negative
result. Likewise an absent HostState component remains unknown. Feature
provenance is carried by its unchanged manifest and output location; no feature
values or evidence are synthesized. A genome-derived annotation may have been
computed after HostState `as_of`, but must exist before bundle creation. It does
not mutate H(t) or establish a contemporaneous dynamic state.

`content_sha256` hashes canonical JSON (UTF-8, sorted keys, compact separators),
excluding only the outer hash field. Embedded snapshots and resource byte hashes
are covered. PhageProfile is hashed in full. Sequence normalized hashes use SDL's
uppercase DNA contigs, ignoring headers/whitespace, with a NUL before each contig;
raw file hashes cover exact bytes including headers. These are different digests.
Existing export directories cannot be overwritten. New inputs require a new
bundle identity/export. Relocating the whole directory preserves validity.

## Results

PairPrediction echoes exact `bundle_id` and `bundle_sha256`. The endpoint, pair,
context and all inputs are thereby bound together. Results cannot predate the
bundle. Model output is one of:

- `not_estimated`: null value and uncertainty, with an explicit reason.
- `uncalibrated_score`: a bounded score; model and preprocessing artifact digests
  are mandatory. It is not a calibrated infection probability.
- `probability`: additionally requires a verified source assay definition,
  calibration artifact digest, containing uncertainty interval and named method.

These checks establish consistency, not scientific verification of a calibration
claim. Model promotion requires the separate evidence and evaluation review.
Unresolved outputs never become training negatives.

## Runnable round trip

With SDL and model checkouts installed in the same Python environment, start in
SDL. The checked-in bundle is a tiny software-only fixture, not experimental data:

```bash
python -m phage_lab validate-model-bundle tests/fixtures/input_bundle/bundle.json --require-sequence host
phage-model-interface validate-bundle tests/fixtures/input_bundle/bundle.json --require-sequence host
phage-model-interface create-unresolved-pair-prediction tests/fixtures/input_bundle/bundle.json /tmp/pair-unresolved.json --prediction-id interface-check
python -m phage_lab validate-pair-prediction /tmp/pair-unresolved.json --bundle tests/fixtures/input_bundle/bundle.json
```

The output path must be new and outside the frozen bundle directory. This round
trip runs no inference. A model implementation uses `load_input_bundle`, performs
its own preprocessing/inference, and calls `build_pair_prediction` with explicit
outputs and version digests. `read_sequence_bytes` returns verified FASTA bytes.

For real inputs SDL's `export-model-bundle` takes `--host-state`, `--phage-profile`,
`--endpoint`, `--assay-context`, `--bundle-id` and `--outdir`. Optional
`--sequence-registry` plus repeated `--include-sequence host|phage` explicitly
copies sequences. Repeated `--feature-artifact` includes validated annotation
outputs. See SDL `docs/usage.md` for an export example.

Cross-check both repositories from SDL:

```bash
python tests/integration/check_model_bundle.py --model-repo ../phage-host-interaction-predictor
```

## Mechanistic targets

The earlier ModelRequest v1 / ModelResult v1 interface (adsorption, productive
infection and resistance emergence, bound by identity only) has been removed from
this repository; it carried no usable input data.

**A source-interaction PairPrediction is not automatically converted to a
mechanistic target or admitted to SDL's productive-infection CandidateSet policy.**
That needs a reviewed endpoint mapping and an endpoint-specific policy. The first
label remains unverified, so this revision deliberately provides no such mapping.
Do not copy a source score into a mechanistic probability slot.

## Training and feedback

SDL owns QC-admitted physical results and explicit assessor labels. Curated
training-data export is separate from inference-input export. Training pipelines
must exclude target outcomes and future evidence from features, group duplicate
isolates/genomes, and keep preprocessing/calibration off the held-out evaluation
partition. A frozen input hash alone does not prove freedom from label leakage.
Retraining and promotion are explicit model-side actions, never automatic
consequences of a new laboratory result.
