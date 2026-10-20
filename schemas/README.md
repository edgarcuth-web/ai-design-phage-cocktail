# Model exchange schemas

`model-request-v1.schema.json` and `model-result-v1.schema.json` define the
versioned wire boundary between this repository and
`self-driving-lab-for-phage-therapy`. The same schema bytes are mirrored in the
SDL repository and checked by tests.

Scientific characterization, sequence registries, FeatureArtifact producers and
PADLOC implementations are owned by SDL. The bundle schema includes the exported
FeatureArtifact wire shape so a consumer can validate it without importing SDL. Model-specific representations, preprocessing
configurations, training/evaluation records, checkpoints and calibration artifacts
belong here.


`model-input-bundle-v1.schema.json` and `pair-prediction-v1.schema.json` define
the current source-endpoint exchange. Bundle resources carry readable payloads;
legacy request hashes alone do not. These schemas and
`tests/fixtures/input_bundle/` are mirrored in both repositories. Embedded
FeatureArtifact definitions mirror SDL's owned v1 schema; no analyzer code moves
into the model repository. The cross-repository check is documented in
`docs/model-interface.md`.
