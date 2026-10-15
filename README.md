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

## First milestones

- [x] Create a reproducible Python project skeleton
- [x] Implement a naive independence-based cocktail baseline
- [ ] Import a small public phage-host interaction dataset
- [ ] Explore class balance and train/test splitting
- [ ] Train a single-phage baseline
- [ ] Calibrate pairwise infection probabilities
- [ ] Evaluate cocktails on unseen bacterial strains
- [ ] Compare against stronger phage-host prediction baselines

## Scope

This repository is a research prototype for computational prediction. It does not establish clinical efficacy or substitute for experimental validation.
