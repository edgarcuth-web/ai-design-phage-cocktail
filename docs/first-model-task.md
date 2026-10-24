# First model task: phage–host source-interaction prediction

Status: task, exchange and a fixed host-feature development model are implemented.
The first study did not meet the criterion for improving over phage prevalence;
see [its audit and results](host-baseline-results.md). Calibrated inference has not
been demonstrated. Project name: **Phage–Host Interaction Predictor**. The existing
Python package name remains `phage_cocktail` so current scripts keep working.

The later [paired baseline study](paired-baseline-results.md) adds explicit phage
4-mer representations and frozen cold-phage/double-cold evaluation. Its primary
improvement criterion was not met; it does not redefine this first task's endpoint.

## Task and evidence

Given one bacterial strain's available representation and one existing phage's
identity or representation, predict the binary label of a named, pinned source
dataset. Batch calls form a pair matrix. They do not estimate mixed-population
effects, cocktail efficacy, or the properties of a hypothetical engineered phage.

| Item | First version |
| --- | --- |
| Inputs | Exact HostState and PhageProfile snapshots; supplied scientific FeatureArtifacts; optional, explicitly exported FASTA resources with byte and normalized-sequence hashes |
| Endpoint | `genophi-ecoli-subset:source-interaction:v1`, source column `interaction` |
| Labels | Source-coded 1/0; `infection` is only a rename. The later primary-source audit excludes 41 zeros with missing MLC; absent/missing observations remain unknown |
| Assay context | Explicitly unknown for the current subset; caller-declared context must include evidence |
| Output | A source-endpoint score, or `not_estimated` with a reason; a probability additionally requires a verified endpoint, calibration artifact and uncertainty |
| Generalization | Unseen bacterial strain IDs with phages that may be known in training; not unseen phages, species, laboratories or conditions |
| Primary comparator | Training-only, Beta(1,1)-smoothed per-phage prevalence; also report global prevalence and 0.5 baselines |

The endpoint descriptor in
[`data/endpoints/`](../data/endpoints/genophi_ecoli_source_interaction.v1.json)
is derived from the dataset registry. Changing its definition requires a new
reviewed endpoint/version and creates different bundle/result hashes. It does
not map to the legacy `productive_infection` target.

FASTA input means already admitted, assembled DNA. FASTQ QC/assembly and scientific
annotation are SDL ingestion tasks and are not implemented by this model
interface. Genome models must request the sequence inputs they actually require;
`load_input_bundle(..., required_sequences=("host", "phage"))` rejects missing
payloads. Models using existing phage IDs must state that dependence in their
model card and abstain on unsupported phages.

## Evidence review, 2026-10-19

| Evidence | What it establishes | Remaining gap |
| --- | --- | --- |
| [Pinned GenoPHI matrix](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/data/interaction_matrices/ecoli_interaction_matrix_subset.csv) and local registry | Exact downloadable bytes and expected 177 × 94 = 16,638 binary entries | Entries alone do not establish primary observation status |
| [Pinned upstream README](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/README.md) | Calls this an E. coli subset for testing; describes binary phenotype input | Does not specify this subset's assay threshold or recoding rule |
| [Noonan et al., 2026](https://doi.org/10.1038/s41564-026-02482-5), model-development datasets | Identifies the same subset dimensions and cites its source | Its separate new-phage validation assay cannot establish the old subset's label mapping |
| [Gaborieau et al., 2024](https://doi.org/10.1038/s41564-024-01832-5) and its pinned source matrix | All 16,597 known readings match MLC > 0; 41 unknowns were source-coded zero | Raw-image/replicate aggregation and row-specific assay contexts were not revalidated |

The original registry/endpoint remains frozen and unverified; the new versioned
study explicitly records the primary-source reconciliation and exclusions.
Mechanistic target mapping remains unverified. The upstream software's MIT
statement does not independently resolve dataset reuse terms.

To close the label review, retain the source table/protocol reference, original
values, exact binary transformation, ambiguous/missing exclusions, subset ID
mapping and reproducible row-level agreement report. Do not infer these from a
publication title or its use of the word “infection”.

## Acceptance

Software acceptance is executable now:

1. Both repositories accept identical bundle and result fixtures; every supplied
   file is readable and hash/size checked. Relocation preserves validity.
2. Changed bytes, wrong pair identity, missing required inputs, path escapes,
   future inputs and wrong-bundle results fail. Partial/absent annotation remains
   partial/absent; no zeros are synthesized.
3. Scores cannot silently become calibrated probabilities or legacy mechanistic
   predictions. No inference produces an explicit `not_estimated` result.

Scientific acceptance for a **host-aware model** is separate:

1. Verify the registered data checksum, source endpoint and host/phage ID mapping.
   Freeze dataset, representation, preprocessing, code, weights and split IDs.
2. Preserve the historical strain-ID benchmark for reproduction. New training
   uses observed-label curation and genomic host groups (70/15/15 groups, seed 42).
   Fit preprocessing on training data; fit thresholds/calibration on validation
   only. A strain-ID split alone does not prove genomic separation.
3. Report AUROC, average precision, Brier and log loss against all three baselines
   on the same observed rows, with counts and missingness. The primary gain is
   average precision over per-phage prevalence, with paired host-group bootstrap
   95% intervals. Require its lower bound to exceed zero before claiming improved
   discrimination; otherwise report the inconclusive/negative finding.
4. Probability promotion additionally needs verified assay semantics, documented
   calibration/uncertainty, no worse Brier/log loss than the comparator on the
   locked evaluation cohort, and a model-card review. Hash admission is not that
   review. Reserve a new independent holdout for release because this development
   test set's results have already been inspected.
5. Evaluate unseen phages and experimental-domain shifts separately. Do not claim
   those capabilities from the first split. Physical feedback is a new curated
   dataset version, never an automatic model update.

## Ownership

SDL owns sequence ingestion, shared biological annotation, evidence/state,
retrieval, candidate decisions, experiments and feedback. This repository owns
model-specific transformations, training, evaluation, calibration and inference.
The first delivery module is **Candidate Phage Prioritization** in SDL; the
predictor is one dependency. Existing source-label scores remain research outputs
until an endpoint-specific selection policy is validated. The legacy SDL
CandidateSet policy consumes productive-infection estimates and cannot consume
these source-label scores by relabelling them.

The old independence/coverage functions remain offline research baselines. Their
union of individually observed positive pairs is not a measured cocktail outcome.
