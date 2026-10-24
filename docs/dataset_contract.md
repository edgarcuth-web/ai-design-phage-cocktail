# Dataset and label contract v1

## Registry is the source of truth

Each dataset used for training or evaluation must have a versioned JSON manifest
under `data/registry/`. The manifest records the upstream source revision,
download URL, local path, SHA-256, schema mapping, expected summary statistics,
license-review status, and label semantics.

A dataset is usable only if its bytes match the registered SHA-256 and its
validated summary matches the manifest. Updating upstream data requires a new
manifest version or an explicit reviewed manifest change; do not silently
replace a file behind an existing dataset ID.

## Current dataset: GenoPHI E. coli subset

The source uses columns `strain, phage, interaction`. The canonical table uses
`bacterium, phage, infection` for compatibility with existing code. The rename
is syntactic only.

The current manifest intentionally records the operational assay definition and
mapping to mechanistic model targets as **not verified**. In particular, a
binary source interaction label must not automatically be treated as a direct
measurement of any one of:

- adsorption;
- productive infection;
- resistance emergence.

Those are distinct targets in the lab/model exchange contract. A future target
mapping requires reviewing the source phenotype definition and assay protocol,
then documenting the rationale and uncertainty.

## Missingness and negatives

An explicit primary MLC reading of `0` is an observed source negative. An absent
or missing primary reading is **unknown/unobserved**, not a negative. A binary CSV
entry alone does not prove a primary reading existed: the row-level audit found
41 source zeros corresponding to missing MLC. New training uses the study's
curated observed rows; the historical registry preserves the original bytes.
See [the audit and results](host-baseline-results.md).

Pairwise metrics are computed only on explicitly observed test rows.

## Evaluation split

The host-feature study assigns whole genomic host groups, defined by published
core distances and canonical assembly content, 70%/15%/15% by group with seed 42.
Its protocol and exact split IDs are frozen under `data/studies/`; group fractions
need not equal host fractions. All baselines are recomputed on those same observed
evaluation rows. Phage identities may be known during training, so this split
does not test unseen phages; the [paired split freeze](paired-split-freeze.md)
adds frozen cold-phage and double-cold regimes.

## Reproducibility

Before training or evaluation:

1. resolve the dataset manifest;
2. verify the raw-file SHA-256;
3. validate the canonical schema and duplicate-pair constraints;
4. verify expected summary statistics;
5. reconcile primary label missingness and freeze deterministic host-group splits;
6. record dataset ID/version, manifest, split seed, code revision, model artifact,
   and evaluation metrics together.
