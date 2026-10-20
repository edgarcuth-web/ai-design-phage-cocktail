# Host baseline study v1: frozen development protocol

The executable protocol is
[`data/studies/genophi_host_baseline.v1.json`](../data/studies/genophi_host_baseline.v1.json).
It is written before fitting or examining this study's model metrics. Its hash is
recorded with every run. Model selection is not performed against validation or
test results in this study.

The registered 177-host, 94-phage subset is an index of candidate pairs. Published
MLC observations, not the subset's filled binary matrix, determine which pairs
are admissible. Explicit MLC zero is negative; scores 1–4 are positive; missing
MLC remains unknown and is excluded with a row-level reason. The old registered
dataset and its results remain a separate development record.

The model uses supplied biological annotations: DefenseFinder subtype counts and
four surface/capsule categories. It does not run annotation tools. In particular,
the published metadata's `n_infections` column is an outcome-derived quantity and
is never an input. Phage identity chooses a separately fitted logistic model;
this study cannot extrapolate to an unseen phage.

Training-only host preprocessing and fixed L2 logistic regression (C=0.1) are
compared with 0.5, global prevalence and per-phage Beta(1,1) prevalence. Genomic
near-neighbours (published core distance <0.0001) and identical assembly content
are grouped before deterministic 70/15/15 group assignment. The threshold comes
from the original study's published grouping rationale; the supplied distance
matrix is used only for partitioning, never as a fitted feature.

The primary contrast is test average precision minus per-phage prevalence on
exactly the same observed rows. A paired bootstrap resamples held-out host groups,
keeping their full observed phage panels together. It provides a conditional
development interval, not uncertainty across laboratories or new phages. Report
AUROC, Brier and log loss too; scores are not promoted to calibrated probabilities.
The already inspected development data cannot supply an independent release
evaluation.
