# Model exchange schemas

The canonical v1 wire schemas are `model-request-v1.schema.json` and `model-result-v1.schema.json`. The lab repository mirrors these files and uses the same conformance fixtures. Update both copies together; cross-record hashes, lineage, timestamps, and target-scope admission remain enforced by the lab's semantic contract.


The model repository also defines the internal `feature-artifact-v1.schema.json` provenance manifest for computational feature outputs. It is not a ModelRequest/ModelResult wire schema and does not imply that a feature extractor is installed or has run. Its packaged copy under `src/phage_cocktail/schemas/` must remain byte-identical.
