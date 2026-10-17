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

An explicit source label of `0` is an observed negative label according to the
source's coding. An absent bacterium-phage row is **unknown/unobserved**, not a
negative. Loaders and evaluation code must not impute absent rows as zero.

Pairwise metrics are computed only on explicitly observed test rows. Whole-cohort
cocktail coverage requires a complete observed test matrix for every selected
phage and every held-out strain; otherwise the coverage calculation fails closed
instead of shrinking the denominator silently.

## Evaluation split

The default split is 70%/15%/15% by bacterial strain, deterministically assigned
with seed 42. A strain belongs to exactly one of train, validation, or test.
Benchmark functions independently reject train/test strain overlap.

This split tests prediction on **unseen bacterial strains with phage identities
that may be known during training**. It does not test cold-start generalization to
unseen phages, nor does it establish clinical performance. Later evaluations
should add group-aware repeated splits, uncertainty intervals, and a separate
unseen-phage protocol where scientifically appropriate.

## Reproducibility

Before training or evaluation:

1. resolve the dataset manifest;
2. verify the raw-file SHA-256;
3. validate the canonical schema and duplicate-pair constraints;
4. verify expected summary statistics;
5. create the deterministic strain split;
6. record dataset ID/version, manifest, split seed, code revision, model artifact,
   and evaluation metrics together.
