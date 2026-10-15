# Phage Cocktail Predictor

Predict whether a bacteriophage or phage cocktail can successfully infect a bacterial strain, while making uncertainty explicit when the target bacterium is incompletely characterized.

## MVP question

Given a bacterial strain and a set of candidate phages:

1. estimate pairwise infection probabilities for each phage-host pair;
2. combine those probabilities into a simple cocktail baseline;
3. evaluate on bacterial strains not seen during training;
4. report calibrated uncertainty rather than only a binary prediction.

The first cocktail baseline assumes independent phage action:

```text
P(at least one phage infects) = 1 - product(1 - p_i)
```

This is deliberately simple. Later versions can learn synergy, antagonism, shared receptors, bacterial defence systems, and other dependencies.

## Repository layout

```text
data/
  raw/            # untouched source data
  processed/      # cleaned/model-ready data
notebooks/        # exploration and experiments
src/phage_cocktail/
  baseline.py     # simple cocktail probability baseline
  data.py         # interaction-table loading/validation
tests/            # automated tests
```

## Expected interaction-table format

```csv
bacterium,phage,infection
Ecoli_01,P001,1
Ecoli_01,P002,0
Ecoli_02,P001,0
```

`infection` must be 0 or 1.

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

## Step 1: fetch and inspect the first dataset

Download the pinned GenoPHI E. coli subset:

~~~bash
python scripts/fetch_genophi.py
python scripts/explore_interactions.py
~~~

For interactive exploration, open:

~~~text
notebooks/01_explore_interactions.ipynb
~~~

## Step 2: leakage-safe benchmark

Run the fixed held-out-strain benchmark:

~~~bash
python scripts/run_benchmarks.py
~~~

Or open:

~~~text
notebooks/02_leakage_safe_benchmark.ipynb
~~~

The split is 70%/15%/15% by bacterial strain, so no target strain is shared
between train, validation, and test. Benchmark results and interpretation are
recorded in docs/step2_benchmark.md.

## Step 3A: reproduce a genome-aware reference

Before building our own genome model, reproduce GenoPHI's official 25-strain x
25-phage protein-family demo. This verifies the complete path from protein
FASTA files through MMseqs2-derived features to strain-separated prediction.

On Linux/macOS, or inside WSL2 on Windows:

~~~bash
conda env create -f environment.yml
conda activate phage-cocktail

python scripts/bootstrap_genophi_demo.py
python scripts/run_genophi_reference.py
python scripts/summarize_genophi_reference.py
~~~

GenoPHI is pinned to version 1.0.3 and its demo data are pinned to commit
f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f. Details are in
docs/step3_genophi_reference.md.

The pinned matrix contains 16,638 interactions across 177 bacterial strains and
94 phages, with a 22.9% positive interaction rate. See data/README.md for
provenance.

## First milestones

- [x] Create a reproducible Python project skeleton
- [x] Implement a naive independence-based cocktail baseline
- [x] Import a reproducible public phage-host interaction dataset
- [x] Explore class balance and matrix coverage
- [x] Define leakage-safe held-out-strain train/validation/test splitting
- [x] Establish host-agnostic pairwise and broadest-cocktail baselines
- [x] Wire in the official GenoPHI genome-aware reference workflow
- [ ] Reproduce GenoPHI demo metrics locally
- [ ] Train our own host-aware single-phage model on the fixed benchmark
- [ ] Calibrate pairwise infection probabilities
- [ ] Evaluate cocktails on unseen bacterial strains
- [ ] Compare against stronger phage-host prediction baselines

## Scope

This repository is a research prototype for computational prediction. It does not establish clinical efficacy or substitute for experimental validation.
