# Model exchange schemas

`model-request-v1.schema.json` and `model-result-v1.schema.json` define the
versioned wire boundary between this repository and
`self-driving-lab-for-phage-therapy`. The same schema bytes are mirrored in the
SDL repository and checked by tests.

Scientific characterization and its schemas—including sequence registries,
FeatureArtifact provenance, and PADLOC result summaries—are owned by SDL and are
not part of this model repository. Model-specific representations, preprocessing
configurations, training/evaluation records, checkpoints and calibration artifacts
belong here.
