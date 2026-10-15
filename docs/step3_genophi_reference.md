# Step 3A: genome-based reference model

Before implementing a custom genome model, reproduce the official GenoPHI
protein-family workflow on the authors' small end-to-end demo.

This is an external reference baseline, not the final model in this repository.

## Why this step exists

Step 2 showed that phage identity alone is already a strong nuisance baseline.
A genome-aware model must therefore show incremental predictive value on
unseen bacterial strains.

GenoPHI is a useful reference because it:

- derives features from bacterial and phage protein sequences;
- clusters proteins into protein families with MMseqs2;
- performs feature selection;
- trains CatBoost classifiers;
- supports strain-separated evaluation with filter_type=strain.

The demo contains 25 E. coli strains, 25 phages, and 625 measured
strain-phage interactions.

## Pinned versions

GenoPHI package:

    1.0.3

GenoPHI source/data commit:

    f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f

The repository checkout is sparse: only data/test_data is required.

## Environment

MMseqs2 is an external dependency. The upstream project documents Linux and
macOS as tested operating systems. On Windows, use WSL2/Linux for this
reference workflow.

With conda or mamba:

    conda env create -f environment.yml
    conda activate phage-cocktail

Or install the Python portion into an existing environment:

    pip install -e ".[dev,reference]"

and install MMseqs2 separately.

Verify:

    genophi --version
    mmseqs version

## Fetch the pinned demo data

    python scripts/bootstrap_genophi_demo.py

This creates:

    external/GenoPHI/data/test_data/
        strain_AAs/
        phage_AAs/
        ecoli_test_interaction_matrix.csv

The external directory is ignored by Git.

## Run the authors' demo settings

    python scripts/run_genophi_reference.py

The command uses:

    --num_features 50
    --num_runs_fs 5
    --num_runs_modeling 10
    --method rfe
    --filter_type strain

Results are written under:

    artifacts/genophi_demo/

The artifacts directory is ignored by Git.

## Summarize the result

    python scripts/summarize_genophi_reference.py

The summarizer looks for GenoPHI's model-performance CSV and reports mean and
standard deviation for classification metrics such as AUC, accuracy,
precision, recall, F1 and MCC when present.

## Important limitation

The 25 x 25 demo is for verifying the complete genome-to-prediction pipeline.
It is too small to replace the 177 x 94 benchmark from Steps 1-2.

After the reference workflow runs successfully, the next task is to connect
genome-derived features to our fixed held-out-strain benchmark and compare
against the Step 2 phage-prevalence baseline on exactly the same split.
