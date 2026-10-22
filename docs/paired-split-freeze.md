# Paired benchmark: source review and split freeze

The [v2 plan](../data/studies/paired_benchmark.v2.json) builds on the unchanged
v1 proposal and preflight. This delivery freezes operational phage-neighbour
groups and all nine development split manifests. No model has been fitted,
scored, calibrated or admitted to SDL.

## Source review, 2026-10-23

Two additional files were pinned to the same original-study revision:
the [raw annotations](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/data/interactions/raw/raw_interactions.csv)
and [phage metadata](https://github.com/mdmparis/coli_phage_interactions_2023/blob/99656d8c020fcd1dbbc30a02ee10274a07655562/data/genomics/phages/guelin_collection.csv).
Their exact byte sizes, SHA-256 and native Git blob identities are in the
[source lock](../data/studies/paired_source_review.resources.v1.json).

All 94 locked phage genome files match their native source Git blobs and exact
record headers. The five duplicate groups are already present under different
names in those source files. All five also have the same recorded species and
isolation host within each pair. This supports a correct source-file join;
it does not establish identical physical stocks or exclude an upstream error.

The [source-review report](../data/studies/results/paired_split_freeze.v1/source_review.json)
traces the 65 conflicting sequence pairs into 318,816 original annotation rows:

| Trace result | Count |
| --- | ---: |
| Original annotation rows linked to the conflicting records | 1,134 |
| Matched nominal contexts (host, image, replicate, plate, dilution) | 567 |
| Matched contexts with interpretable scores on both records | 551 |
| Different interpretable scores in a matched nominal context | 147 |
| Sequence pairs with at least one such difference | 61 of 65 |
| Sequence pairs missing rows in this raw table | 2 |
| Remaining pairs without a directly reproduced interpretable contrast | 2 |

The missing records are `536` and `LF82`, each paired with `BCH953_P4` and
`BCH953_P5`. Both hosts appear as isolation strains in source metadata, consistent
with the original README's description of separate assessment for isolation
strains. This is a source-coverage explanation, not verification of those missing
measurements. Raw `n` remains uninterpretable, never a negative.

Images were not independently rescored. The mapping from individual annotations
and replicates to aggregate MLC was not rebuilt. Matching nominal metadata does
not establish matching physical concentration, stock history or host state. No
confirmed cause, corrected target or majority-vote consensus is created. All
16,597 observed aggregate source labels and 41 missing exclusions remain intact.
Source images, replicate scores, MLC and audit diagnostics are excluded as model
input features.

## Operational sequence-neighbour grouping

The [grouping definition](../data/studies/phage_grouping.v1.json) fixes
minimap2 2.31-r1302, source commit
`3c28777e7e2dcc90f825de1b9f17a89cca7d4452`, and the complete command options.
The run records the actual executable SHA-256. The official versioned
[manual](https://github.com/lh3/minimap2/blob/v2.31/minimap2.1) defines PAF counts,
base-level CIGAR output and all-chain mapping. This study uses `asm20`, `-c`,
`-P`, seed 42 and two threads on the same 94 locked single-contig assemblies as
both references and queries. A different contig scope needs a new implementation.

Our grouping rule is an evaluation choice, not a published taxonomic cutoff:

1. Always connect exact canonical assembly duplicates.
2. For each ordered genome pair, retain base-aligned blocks above the identity
   threshold. Sort by decreasing matches/block length and deterministic coordinates.
3. Greedily retain blocks with no overlap on either genome. Repeated or overlapping
   hits cannot inflate coverage; partial matches need not be collinear.
4. Connect the pair when either mapping direction covers at least 85% of **both**
   complete genomes. Use the connected components, including transitive links.

| Identity threshold | Purpose | Components | Undirected edges |
| --- | --- | ---: | ---: |
| 95% | Primary | 69 | 30 |
| 90% | Prespecified broader-neighbour sensitivity | 50 | 91 |
| 99% | Prespecified tighter-neighbour sensitivity | 84 | 10 |

Identity is matches divided by alignment block length including gaps. The rule
does not guarantee that all members of one component meet the threshold pairwise.
A missing heuristic alignment edge does not prove unrelatedness. These groups
do not exclude shared proteins/RBP domains, remote ancestry or pretraining
overlap, and are not interpreted as species/genus labels or functional identity.
No infection label or model metric selects thresholds, edges or random seeds.

## Frozen rows and feasibility

The [freeze report](../data/studies/results/paired_split_freeze.v1/report.json)
binds the grouping definition, input preflight, source locks, executable, command,
implementation and all generated artifact hashes. Committed compressed PAF
alignments, group/edge CSVs and row-ID manifests support independent inspection.
Gzip output has a fixed timestamp; input row order does not affect membership.
No genome sequence payloads or row-level target labels are committed here.

Hosts retain the 147-group rule. Host and phage groups are assigned independently;
fit uses H_train × P_train minus the sequence-pair-hashed warm masks. All observed
rows have exactly one fit, validation, test or unused assignment. Identical sequence
pairs share their warm mask. Intended seen entities must remain represented in
fit, and fit must retain both source-label classes.

The primary scheme is 95% identity, bilateral coverage 85%, seed 42:

| Cell | Rows | Positive labels | Host groups | Phage groups |
| --- | ---: | ---: | ---: | ---: |
| Fit | 6,834 | 1,493 | 103 | 48 |
| Warm-pair validation | 374 | 77 | 98 | 48 |
| Warm-pair test | 370 | 81 | 95 | 48 |
| Cold-host validation | 1,690 | 305 | 22 | 48 |
| Cold-host test | 2,210 | 572 | 22 | 48 |
| Cold-phage validation | 1,746 | 508 | 103 | 11 |
| Cold-phage test | 1,633 | 314 | 103 | 10 |
| Double-cold validation | 390 | 104 | 22 | 11 |
| Double-cold test | 476 | 108 | 22 | 10 |
| Unused crossed validation/test rectangles | 874 | 251 | 44 | 21 |

All three thresholds and seeds 42/43/44 were frozen. They are not alternative
test sets to choose after inspecting results. Counts establish feasibility,
not model performance; test outcomes remain part of an already inspected
development source, not independent release evidence.

## Reproduce

After running the existing input preflight, supply the pinned minimap2 executable:

```bash
python scripts/freeze_paired_splits.py \
  --preflight artifacts/paired_preflight_v1 \
  --data-root data/raw/genophi_host_study \
  --minimap2 /path/to/minimap2 \
  --fetch --outdir artifacts/paired_split_freeze_v1
```

`--fetch` retrieves only missing locked resources; omit it for a cached offline
run. Existing output directories are refused. All grouping/split artifacts are
derived from admitted inputs, with no model fitting or biological label upgrade.
The original v1 plan/report remain unchanged. The next model implementation must
consume these frozen rows and the v2 selection/weighting rules. Feature extraction
and paired fitting are not implemented by this command; dataset reuse terms remain
a separate documented requirement, not resolved by a source-file join.
