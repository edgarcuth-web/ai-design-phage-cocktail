# Phage–Host Interaction Predictor

Research infrastructure for phage-host interaction prediction and phage-cocktail
baselines. Dataset provenance, label semantics and evaluation boundaries are
explicit; this repository does not establish clinical efficacy.

## First task

Predict a named source-defined binary outcome for one bacterial strain and one
existing phage, then evaluate on unseen strains with possibly known phages.
The task, label evidence, generalization limits and acceptance criteria are in
[`docs/first-model-task.md`](docs/first-model-task.md).

Current delivery: a working file-based model boundary, a real-data label/genome
audit, frozen host-feature weights and a completed additive host–phage benchmark.
The first host model did not establish improvement over phage prevalence. The new
4-mer pair baseline also did not meet its cold-phage improvement criterion across
the prespecified sensitivities. Both remain uncalibrated research scores; see the
[host study](docs/host-baseline-results.md) and
[paired study](docs/paired-baseline-results.md).
Pairwise matrices and offline union coverage do not establish cocktail efficacy.

| Capability | Current status |
| --- | --- |
| File exchange | Bundle/resource validation, result validation and explicit abstention |
| Trained host baseline | Host features plus known phage identity; uncalibrated source-MLC scores; no demonstrated improvement |
| Baseline serving through PairPrediction | Not implemented; weights and the generic exchange are separate paths |
| Mechanistic predictions / unseen phages | Not established by the current baseline |
| Pair representations and pretrained genomic models | Future benchmark candidates; no backbone selected |

## Repository layout

```text
data/
  raw/              # downloaded source data; not committed
  processed/        # cleaned/model-ready data; not committed
  registry/         # versioned dataset manifests and label semantics
  endpoints/        # explicit source-label endpoint descriptors
  studies/          # frozen protocols, resource locks and audit reports
models/             # frozen development JSON weights and model card
notebooks/          # exploration and experiments
src/phage_cocktail/
  artifacts.py      # canonical JSON hashes and raw-file digests
  model_interface.py # exchange validation and CLI
  input_bundle.py   # frozen input/resource consumer and PairPrediction
  baseline.py       # simple cocktail probability baseline
  data.py           # observed interaction-table loading/validation
  registry.py       # dataset checksum and summary verification
  split.py          # deterministic strain-level splits
  benchmark.py      # leakage-guarded pairwise/cocktail baselines
  cohort.py         # published MLC reconciliation and genomic host grouping
  host_baseline.py  # train-only preprocessing and known-phage logistic baseline
  host_study.py     # locked real-data audit/training/evaluation runner
tests/              # automated tests
schemas/            # bundle/pair-result and legacy request/result wire schemas
```

## Dataset and labels

The first dataset is the pinned GenoPHI E. coli subset. The fetch script verifies
the registered SHA-256, and analysis scripts verify the checksum and summary
statistics before use.

The canonical `infection` column represents a source-coded binary label. A later
row-level audit found **41 zeros whose primary MLC readings are missing**. The
new study excludes these, retaining 16,597 observed labels that exactly match
`MLC > 0`. The target combines source spotting-assay categories; a mapping to
adsorption, productive infection or resistance emergence remains unverified.
Spotting/clearing, plaque formation, productive infection and growth suppression
are not interchangeable labels; [endpoint semantics](docs/endpoint-semantics.md)
defines their review boundaries without changing the frozen source endpoint.
See [the audit](docs/host-baseline-results.md),
[`docs/dataset_contract.md`](docs/dataset_contract.md) and
[`data/README.md`](data/README.md).

## Model exchange

SDL exports **ModelInputBundle v1** with complete pair snapshots, an explicit
source endpoint/context, scientific FeatureArtifact provenance and optional local
FASTA payloads. The consumer verifies all supplied bytes and rejects missing
model-required sequences. **PairPrediction v1** binds the exact bundle and
separates unresolved, uncalibrated-score and probability outputs.

```bash
phage-model-interface validate-bundle tests/fixtures/input_bundle/bundle.json --require-sequence host
```

The fixture is synthetic and tests software only. The legacy ModelRequest /
ModelResult v1 commands remain for existing consumers; their three mechanistic
targets are not aliases for source labels. See
[`docs/model-interface.md`](docs/model-interface.md) for the complete round trip.

## Setup

Requires Python 3.11+.

```bash
python -m venv .venv
```

Activate the environment, then install:

```bash
pip install -e ".[dev]"
pytest
```

## Step 1: fetch and inspect the registered dataset

```bash
python scripts/fetch_genophi.py
python scripts/explore_interactions.py
```

For interactive exploration, open `notebooks/01_explore_interactions.ipynb`.

## Step 2: historical binary-matrix benchmark

```bash
python scripts/run_benchmarks.py
```

Or open `notebooks/02_leakage_safe_benchmark.ipynb`.

The split is 70%/15%/15% by bacterial strain. No target strain is shared between
train, validation, and test. Phage identities may be shared across partitions,
so this measures generalization to unseen strains with known phages, not unseen
phages. Cocktail coverage evaluation fails if selected-phage/test-strain pairs
are missing rather than treating missing outcomes as negatives.

These historical results include the 41 source-coded unknowns and are retained
for reproducibility. Use the audited, genomically grouped study below for new
observed-label evaluation. Historical results and their interpretation are in
[`docs/step2_benchmark.md`](docs/step2_benchmark.md).

## First observed-label host-feature study

```bash
python scripts/run_host_baseline.py --fetch --outdir artifacts/host_baseline_v1
```

This verifies 177 host genomes, 94 phage genomes and all annotation/source tables
(approximately 920 MB), excludes unknown MLC readings, freezes host groups and
fits the [fixed protocol](docs/host-baseline-protocol.md). Existing run directories
and tampered input files are rejected. Weights, metrics, row predictions and
provenance are saved together. Frozen development weights are in `models/`.

Test AP is 0.5413 versus 0.4952 for per-phage prevalence; the paired 95% AP-gain
interval is [-0.0229, 0.1092]. This is inconclusive, and Brier/log loss are worse.
See [the full report](docs/host-baseline-results.md) for the split and limitations.

## Paired benchmark protocol

The next paired-model study has a separate
[development protocol and executable input preflight](docs/paired-benchmark-protocol.md).
It audits sequence/label identity and duplicate-genome disagreements without
training. [The follow-up source review and split freeze](docs/paired-split-freeze.md)
now fixes 69 operational phage-neighbour groups for the primary scheme and all
nine development split manifests. It traces source annotation differences without
relabeling them. The [completed fixed baseline study](docs/paired-baseline-results.md)
runs host-only, phage-only and additive pair models, plus applicable comparators,
under both weighting policies. Cold-phage gains remain inconclusive. The already
inspected cohort cannot serve as independent release validation.

## Optional GenoPHI genome-aware reference

GenoPHI's official 25-strain x 25-phage protein-family demo can verify the path
from protein FASTA files through MMseqs2-derived features to strain-separated
prediction.

On Linux/macOS, or inside WSL2 on Windows:

```bash
conda env create -f environment.yml
conda activate phage-cocktail
python scripts/bootstrap_genophi_demo.py
python scripts/run_genophi_reference.py
python scripts/summarize_genophi_reference.py
```

GenoPHI is pinned to version 1.0.3 and its demo data to commit
`f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f`. See
[`docs/step3_genophi_reference.md`](docs/step3_genophi_reference.md).

## Milestones

- [x] Create a reproducible Python project skeleton
- [x] Import a pinned public phage-host interaction dataset
- [x] Register source revision, checksum, expected statistics and label caveats
- [x] Define deterministic held-out-strain splitting and basic baselines
- [x] Guard against train/test strain overlap and missing-pair coverage bias
- [x] Wire in the official GenoPHI genome-aware reference workflow
- [x] Define the first pairwise task and preserve unverified assay semantics
- [x] Consume frozen input bundles and return exact source-endpoint results
- [x] Reconcile source binary labels with published MLC readings and exclude unknowns
- [x] Verify all cohort genomes and supplied biological annotations
- [x] Freeze genomic host groups and run a fixed host-feature baseline
- [x] Report paired host-group uncertainty, including an inconclusive result
- [ ] Finish row-level assay-context and separate dataset-license review
- [ ] Reproduce GenoPHI demo metrics locally
- [ ] Define target-specific label mapping before training mechanistic targets
- [ ] Demonstrate improvement and validate calibration before model promotion
- [x] Evaluate unseen-phage regimes with frozen sequence groups and report inconclusive additive-baseline gains
- [ ] Add repeated group evaluation and an independent release holdout
- [ ] Compare explicit host–phage pair features with available pretrained representations under the same protocol

## Scope

This repository is a research prototype for computational prediction. It does
not establish clinical efficacy or substitute for experimental validation.
