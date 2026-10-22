# Paired benchmark v1: development protocol proposal

The machine-readable plan is
[`paired_benchmark.v1.json`](../data/studies/paired_benchmark.v1.json).
The executable delivery in this change is **input preflight, not model fitting**.
This protocol is a proposal: phage relatedness groups and final split manifests
must be reviewed and locked before any new metrics are inspected. No architecture
or foundation-model checkpoint has been selected.

## Question and target

Does a representation of both the host and the phage improve prediction beyond
host-only susceptibility and phage prevalence, especially for unseen phages?

The target remains observed source `MLC > 0`. Explicit MLC zero is negative;
missing MLC is excluded with a reason. This study does not establish productive
infection, cocktail efficacy or clinical utility. The historical endpoint and
host-baseline study are preserved rather than retroactively redefined.

The locked resources contain 177 host genomes, 94 phage genome records, the subset
index, primary MLC table, host core distances and supplied annotations. The expected
cohort is 16,597 observations and 41 exclusions. Every genome and table must pass
the existing SHA-256/size lock before a full preflight succeeds.

## Executable input preflight

With the repository installed, run from its root:

```bash
python scripts/audit_pair_benchmark.py --fetch --outdir artifacts/paired_preflight_v1
```

For already downloaded resources, omit `--fetch`. `--data-root` can select the
existing host-study download directory. The output directory must be new.
Downloads are of locked public resources; source genomes are not committed.

The audit rechecks source-label agreement, host/phage ID coverage, every sequence
digest, the archived host grouping/split, and the required host annotation join.
It writes observed/excluded rows, host identities/splits, exact phage identities
and all conflicting labels for identical host/phage sequence pairs. It preserves
labels and reports conflicts; it neither merges them nor invents missing pairs.
The report binds protocol, code, resource-lock and output hashes.

`ready_for_training=false` is intentional. Exact-duplicate detection cannot certify
the absence of near-neighbour leakage, the protocol is not yet a training release,
and a successful input audit is not model validation. The preflight emits no model
weights, predicted scores, fitted calibration or biological interpretations.

The observed-row CSV retains the archived host-only `host_group`/`split` columns
for exact cohort reproduction. Those columns are not the new paired-regime split;
the new host/phage split manifests remain a before-fit requirement.

## Preflight findings, 2026-10-22

The saved [preflight report](../data/studies/results/paired_benchmark_preflight.v1.json)
records a fresh verification of all 276 locked files: 177 host genomes, 94 phage
genomes and five source tables. The observed-row, excluded-row and host-split CSV
hashes exactly reproduce the archived host study. No model was fitted or scored.

The 94 phage records reduce to 89 exact assembled sequences. Five duplicate
groups were identified:

| Exact sequence group | Record IDs |
| --- | --- |
| 423_P10 | 423_P10, 423_P9 |
| BCH953_P4 | BCH953_P4, BCH953_P5 |
| BDX03_P1 | BDX03_P1, BDX03_P2 |
| LI10_P1 | LI10_P1, LI10_P2 |
| NAN33_P2 | NAN33_P2, NAN33_P4 |

There are 65 identical host/phage sequence pairs with disagreeing binary labels,
comprising 130 original observations. This is an observed input/label ambiguity,
not proof of a biological mechanism or of annotation error. It motivates retaining
source observations and separating genome grouping from physical stock/context
identity. Per-row details are generated locally by the preflight; the committed
report contains counts, group membership and hashes, not redistributed genomes
or primary observation tables.

## Grouping and split design

Hosts reuse the existing connected components of published core distances below
0.0001 or exact assembled-sequence identity. Recompute from admitted data and
compare the archived 147-group development split. This reproduces a development
boundary; it does not create a new independent holdout.

Exact phage duplicates must share a group even when record IDs differ. The source
contains 94 phage names but only 89 exact assembled sequences. Before fitting,
freeze a sequence-relatedness graph with tool/version, coverage/similarity rules,
parameters, exact-duplicate union, memberships, edge hashes and a planned threshold
sensitivity analysis. Choosing that graph requires a sequence audit, not model
performance. Do not call an ID-only or exact-duplicate-only split lineage held out.

Assign whole host and phage groups independently to 70/15/15 partitions using the
existing deterministic ordering, primary seed 42. Seeds 43 and 44 are prespecified
sensitivity runs, not opportunities to choose a favorable split. Inspect counts
and feasibility before fitting; any revision needs a new recorded protocol version.

| Regime | Validation rows | Test rows | Training relationship |
| --- | --- | --- | --- |
| Warm pair | Fixed masked rows in H_train × P_train | Separate fixed mask in the same rectangle | Both entity groups represented in fit |
| Cold host | H_validation × P_train | H_test × P_train | Host groups absent from fit |
| Cold phage | H_train × P_validation | H_train × P_test | Phage groups absent from fit |
| Double cold | H_validation × P_validation | H_test × P_test | Both entity groups absent from fit |

All models fit the same H_train × P_train rows minus the warm masks. Masks use
the exact sequence-pair hash so duplicate names cannot leak an identical pair
between fit and evaluation. Fail if a group intended to be seen disappears from
fit, a required training label class is absent, or group membership crosses an
excluded axis. The two mixed validation/test rectangles remain unused. Export
every admitted row and exclusion with its regime, partition and exact group IDs.

## Minimal comparisons

| Model | Purpose | Applicable regimes |
| --- | --- | --- |
| Train-only global prevalence | Sanity baseline | All |
| Train-only per-phage prevalence | Known-panel comparator | Seen phages only |
| Existing per-phage host logistic, refitted | Preserve the established baseline | Seen phages only |
| Pooled host-only logistic | Host susceptibility without phage information | All |
| Phage-only 4-mer logistic | Phage compositional signal without host information | All |
| Host + phage additive logistic | Incremental paired information | All |

Use the published 137 defense counts and four surface/capsule categories for the
host block; use canonical within-contig DNA 4-mer frequencies for the phage block.
The JSON defines the proposed extraction and logistic settings. Concatenation is
an additive pair baseline, not an explicit mechanistic interaction model. Entity
ID strings, outcome-derived metadata and evaluation labels are excluded as model
features. Grouping distances are partitioning inputs, not fitted features.

Fit entity-level scaling and category vocabulary on unique training entities only.
Select C from 0.01, 0.1 and 1 using cold-phage validation AP; prefer lower C on ties.
The historical seen-phage-only models keep their fixed configurations. Refit every
applicable comparator on the new fit rows: frozen old weights were trained with
all phage IDs and would violate the new cold-phage boundary. Unsupported models
abstain explicitly; never silently substitute an unseen phage's prevalence.

An optional pretrained-representation comparison is a separate admitted experiment
after checkpoint, tokenizer, sequence windows, pooling, compute budget and known
pretraining-overlap limitations are locked. It is not required to run the minimal
comparison. Mechanism blocks are later grouped ablations on the same rows.

## Evaluation and release boundaries

The primary contrast is cold-phage pooled AP: host+phage additive logistic minus
pooled host-only logistic. Report global/phage-only comparators and all other
regimes separately, with denominators, positive counts and abstention coverage.
Also report AUROC, Brier, log loss and macro per-host/per-phage AP with the number
of undefined single-class groups. A score is not promoted to a calibrated
probability merely because it is bounded or was produced by logistic regression.

Use paired group bootstrap with 2,000 draws and the fixed seed in the JSON.
Resample host groups for cold-host evaluation and phage groups for cold-phage
evaluation, conditional on the seen panel. For double-cold and warm-pair rows,
resample both axes independently and multiply row multiplicities. Report omitted
single-class draws. Fewer than five groups on a resampled axis, missing classes
or no valid draws means descriptive results only. A lower primary AP-gain bound
above zero supports a development improvement for this source, not a biological
mechanism or a release claim. Prespecified sensitivity runs are reported together.

Report sequence-pair label discordance explicitly. Identical assembled inputs
with different labels may reflect record assignment, material history, assay
context or other unrecorded differences; the audit cannot determine the cause.
The primary analysis retains equal weight for every source observation; the
prespecified sensitivity gives each exact sequence pair equal total weight in
both fitting and metrics. Report both; neither analysis turns a conflicting pair
into a voted consensus label.

The source has already been examined during baseline development. Repartitioning
it cannot restore independence. New independently reserved data are required for
release, prospective calibration/OOD assessment and cross-laboratory claims.
Dataset reuse terms also require their own review; software licensing is not a
substitute. No step in this protocol automatically updates SDL decisions or models.
