# Phage Cocktail Predictor

Research infrastructure for phage-host interaction prediction and phage-cocktail
baselines. Dataset provenance, label semantics and evaluation boundaries are
explicit; this repository does not establish clinical efficacy.

## MVP question

Given a bacterial strain and a set of candidate phages:

1. estimate pairwise interaction probabilities for explicitly defined targets;
2. combine those probabilities into a simple cocktail baseline;
3. evaluate on bacterial strains not seen during training;
4. report calibrated uncertainty only when supported by validation.

The first cocktail baseline assumes independent phage action:

```text
P(at least one phage infects) = 1 - product(1 - p_i)
```

This is deliberately simple. Real phages can share receptors, encounter the same
bacterial defence systems, or otherwise have dependent outcomes.

## Repository layout

```text
data/
  raw/              # downloaded source data; not committed
  processed/        # cleaned/model-ready data; not committed
  registry/         # versioned dataset manifests and label semantics
notebooks/          # exploration and experiments
src/phage_cocktail/
  baseline.py       # simple cocktail probability baseline
  data.py           # observed interaction-table loading/validation
  registry.py       # dataset checksum and summary verification
  split.py          # deterministic strain-level splits
  benchmark.py      # leakage-guarded pairwise/cocktail baselines
tests/              # automated tests
schemas/            # versioned ModelRequest / ModelResult JSON Schemas
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

## External model exchange contract

The file-based boundary to the self-driving lab is defined by the v1 JSON Schemas in [`schemas/`](schemas/) and documented in [`docs/model-interface.md`](docs/model-interface.md). Shared fixtures are validated in CI. This is currently a wire-contract milestone, not an inference implementation: no biological prediction is emitted by these fixtures.

## File-based model interface

The model exchange CLI validates the shared `ModelRequest v1` and `ModelResult v1` contracts. It does not run inference. After installing the project, run `phage-model-interface --help` for commands, or see [`docs/model-interface.md`](docs/model-interface.md) for a request/result round trip and explicit unresolved response behavior.

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
- [ ] Verify the source assay definition and license terms
- [ ] Reproduce GenoPHI demo metrics locally
- [ ] Define target-specific label mapping before training mechanistic targets
- [ ] Train and calibrate a host-aware model on the fixed benchmark
- [ ] Evaluate generalization to unseen phages separately
- [ ] Add grouped repeated evaluation and uncertainty intervals

## Scope

This repository is a research prototype for computational prediction. It does
not establish clinical efficacy or substitute for experimental validation.
