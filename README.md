# Phage–Host Interaction Predictor

Research infrastructure for phage-host interaction prediction and phage-cocktail
baselines. Dataset provenance, label semantics and evaluation boundaries are
explicit; this repository does not establish clinical efficacy.

## First task

Predict a named source-defined binary outcome for one bacterial strain and one
existing phage, then evaluate on unseen strains with possibly known phages.
The task, label evidence, generalization limits and acceptance criteria are in
[`docs/first-model-task.md`](docs/first-model-task.md).

Current delivery: dataset/baseline infrastructure and a working file-based input
and result boundary. No trained host-aware model or calibrated inference is
included. Pairwise matrices and offline union-coverage baselines do not establish
cocktail efficacy.

## Repository layout

```text
data/
  raw/              # downloaded source data; not committed
  processed/        # cleaned/model-ready data; not committed
  registry/         # versioned dataset manifests and label semantics
  endpoints/        # explicit source-label endpoint descriptors
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
tests/              # automated tests
schemas/            # bundle/pair-result and legacy request/result wire schemas
```

## Dataset and labels

The first dataset is the pinned GenoPHI E. coli subset. The fetch script verifies
the registered SHA-256, and analysis scripts verify the checksum and summary
statistics before use.

The canonical `infection` column currently represents the source-encoded binary
interaction label. Its exact assay definition and mapping to adsorption,
productive infection, or resistance-emergence targets are **not yet verified**.
An absent pair is unknown/unobserved, not a negative label. See
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

## Step 2: leakage-safe baseline benchmark

```bash
python scripts/run_benchmarks.py
```

Or open `notebooks/02_leakage_safe_benchmark.ipynb`.

The split is 70%/15%/15% by bacterial strain. No target strain is shared between
train, validation, and test. Phage identities may be shared across partitions,
so this measures generalization to unseen strains with known phages, not unseen
phages. Cocktail coverage evaluation fails if selected-phage/test-strain pairs
are missing rather than treating missing outcomes as negatives.

Existing benchmark results and their interpretation are recorded in
[`docs/step2_benchmark.md`](docs/step2_benchmark.md).

## Step 3A: reproduce a genome-aware reference

Before building a custom genome model, reproduce GenoPHI's official 25-strain x
25-phage protein-family demo. This verifies the complete path from protein FASTA
files through MMseqs2-derived features to strain-separated prediction.

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
- [ ] Verify the source assay definition and license terms
- [ ] Reproduce GenoPHI demo metrics locally
- [ ] Define target-specific label mapping before training mechanistic targets
- [ ] Train and calibrate a host-aware model on the fixed benchmark
- [ ] Evaluate generalization to unseen phages separately
- [ ] Add grouped repeated evaluation and uncertainty intervals

## Scope

This repository is a research prototype for computational prediction. It does
not establish clinical efficacy or substitute for experimental validation.
