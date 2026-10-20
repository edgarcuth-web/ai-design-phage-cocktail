# Observed-label audit and first host-aware baseline

The fixed host-feature model has **not demonstrated an improvement** over the
per-phage prevalence comparator. Its test average precision is higher, but the
paired interval includes zero, validation average precision is lower, and test
Brier/log loss are worse. These development weights are not a calibrated model
or an approved candidate-selection policy.

## Source reconciliation

The [pinned GenoPHI subset](https://github.com/Noonanav/GenoPHI/blob/f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f/data/interaction_matrices/ecoli_interaction_matrix_subset.csv)
contains 16,638 source-coded binary rows. They were joined to the
[original published MLC matrix](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/data/interactions/interaction_matrix.csv),
using unique uppercase, hyphen-removed host aliases and exact phage IDs.
Every one of the 177 host IDs has a unique match; no fuzzy matching is used.

| Published MLC | Pairs | Training interpretation |
| --- | ---: | --- |
| 0 | 12,784 | Observed source negative |
| 1 | 1,492 | Observed source positive |
| 2 | 1,041 | Observed source positive |
| 3 | 680 | Observed source positive |
| 4 | 600 | Observed source positive |
| Missing | 41 | Unknown; excluded |

All **16,597 known** binary labels agree with `MLC > 0`. All **41 missing** MLC
values were encoded as zero in the GenoPHI subset: 39 pairs for `370D`, and two
for `ECOR53`. The old registered bytes and historical benchmark are retained;
new training uses the explicitly curated observed cohort. The saved exclusion
list and full host split manifest are in
[`data/studies/results/`](../data/studies/results/).

The study authors' [Viral Host Range DB description](https://viralhostrangedb.pasteur.cloud/data-source/)
identifies a spotting assay, with MLC coding from absent clearing through clearing
at progressively lower phage concentrations. Scores 3 and 4 distinguish plaques
from full lawn clearing at the lowest concentration. Thus the binary target
combines several source assay categories. It is not a separately verified
productive-infection endpoint. The
[source README](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/README.md)
also says isolation strains were assessed separately; uniform row-level assay
conditions are not inferred. This audit reconciles the published aggregate
matrix, without re-scoring images or rebuilding all replicate-level MLC calls.

## Resources and partitioning

The resource lock contains SHA-256 and size for **276 files** (919,683,545 bytes):
177 host DNA assemblies, 94 phage DNA assemblies and five source tables. All
were read and verified. Host assemblies come from
[Figshare article 25941691, version 1](https://doi.org/10.6084/m9.figshare.25941691.v1);
the metadata declares CC BY 4.0. Phage sequences and annotations are pinned to
the original repository revision above. The
[publication archive metadata](https://doi.org/10.5281/zenodo.13831957) also declares
CC BY 4.0, but an incomplete archive transfer was rejected: archive byte
equivalence was not verified. This does not clear the separate GenoPHI software
and data reuse questions by treating its software MIT notice as a data license.

Canonical contig-content hashing found **177 distinct host assemblies** and
**89 distinct phage assemblies** among the 94 phage IDs. Headers, contig order
and strand are ignored; different fragmentation or circular rotations are not
equated. The phage duplicates do not constitute a held-out-phage experiment.

Published core-genome distances below 0.0001 are connected into host groups,
along with any exact assembly matches. This yields **147 host groups**. Groups
are assigned once with seed 42; no group crosses partitions. Fractions apply
to group counts, so host and row fractions differ.

| Partition | Groups | Hosts | Observed pairs | Positive pairs |
| --- | ---: | ---: | ---: | ---: |
| Train | 103 | 117 | 10,957 | 2,473 |
| Validation | 22 | 26 | 2,444 | 491 |
| Test | 22 | 34 | 3,196 | 849 |

This new cohort/split differs from the historical 124/27/26 host split. Baselines
are recomputed on the same new rows; comparison with the old AP of 0.5078 would
confound changes in data, grouping and model.

## Fixed experiment

The model receives **137 published defense-subtype counts**, plus O-antigen,
H-antigen, outer-core LPS and ABC-capsule categories. It uses supplied annotations;
the verified genome files establish input identity and grouping, and are not
consumed as sequence embeddings. Outcome-derived `n_infections` and the other
unselected metadata never enter the feature matrix.

Count scaling and categorical vocabularies are fitted on unique training hosts
only. One L2 logistic regression per known phage uses C=0.1 with no tuning.
`BCH953_P3` has only one training class and explicitly falls back to smoothed
training prevalence. Saved JSON weights can be reloaded by
`predict_host_model`; unsupported phage identities are rejected.

The protocol was frozen before fitting:
`b717e0c5a5c788f19be84d441548fa0a89b111065827df9b39783ea6fecb12bf`.
The [protocol](host-baseline-protocol.md), resource lock, code file hashes,
dependency versions and output digests are recorded in the
[machine-readable report](../data/studies/results/genophi_host_baseline.v1.json).

| Partition | Model | AUROC | Average precision | Brier | Log loss |
| --- | --- | ---: | ---: | ---: | ---: |
| Validation | 0.5 | 0.5000 | 0.2009 | 0.2500 | 0.6931 |
| Validation | Global prevalence | 0.5000 | 0.2009 | 0.1612 | 0.5035 |
| Validation | Phage prevalence | 0.8056 | 0.5328 | 0.1270 | 0.4023 |
| Validation | Host logistic | 0.7969 | 0.5162 | 0.1358 | 0.4247 |
| Test | 0.5 | 0.5000 | 0.2656 | 0.2500 | 0.6931 |
| Test | Global prevalence | 0.5000 | 0.2656 | 0.1967 | 0.5833 |
| Test | Phage prevalence | 0.7585 | 0.4952 | 0.1655 | 0.4987 |
| Test | Host logistic | 0.7735 | 0.5413 | 0.1681 | 0.5215 |

The primary test AP gain is **0.04614**, with paired host-group bootstrap
95% interval **[-0.02290, 0.10923]** (2,000/2,000 valid resamples, seed 20261019).
It fails the predeclared improvement criterion. The interval conditions on this
known phage panel and development split; it is not evidence of generalization
to new phages, laboratories or treatment outcomes. No hyperparameters or inputs
were changed in response to these results.

## Reproduce

From an installed checkout:

```bash
pip install -e ".[dev]" -c data/studies/host_baseline_environment.v1.txt
python scripts/run_host_baseline.py --fetch --outdir artifacts/host_baseline_v1
```

The first fetch is approximately 920 MB. Cached inputs are hash-checked; changed
files and an existing output directory are rejected. The run saves audited
pairs, exclusions, host groups, feature rows, train-only encoder/weights,
predictions, metrics and provenance. Raw genomes and the full training tables
are not committed.

An independent rerun after a DataFrame construction cleanup produced identical
bytes for observed/excluded pairs, host splits, feature values, predictions and
metrics. The cleanup did not alter the model protocol or numerical results.
Further model changes require a new study protocol; an independent release
holdout and endpoint-specific review remain outstanding.
