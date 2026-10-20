# Step 2 benchmark

Historical binary-matrix result. The later [MLC audit](host-baseline-results.md)
found 41 source-coded zeros with missing primary readings. This table is retained
for reproduction, and is not the current observed-label development cohort.

This benchmark uses the pinned GenoPHI E. coli interaction matrix and a
deterministic bacterial-strain split with seed 42. The dataset checksum and
expected summary are verified against the registry before the benchmark runs.

## Split

| Partition | Strains | Interactions | Positive rate |
| --- | ---: | ---: | ---: |
| Train | 124 | 11,656 | 23.14% |
| Validation | 27 | 2,538 | 18.72% |
| Test | 26 | 2,444 | 26.23% |

No bacterial strain occurs in more than one partition. The split holds out
strains, not phage identities: phages may occur in both train and test. Results
therefore measure generalization to unseen strains for phages that may be known
during training, not cold-start generalization to unseen phages.

## Pairwise baselines on test strains

| Baseline | AUROC | Average precision | Brier | Log loss |
| --- | ---: | ---: | ---: | ---: |
| coin_flip_0.5 | 0.5000 | 0.2623 | 0.2500 | 0.6931 |
| global_prevalence | 0.5000 | 0.2623 | 0.1944 | 0.5780 |
| phage_prevalence | 0.7678 | 0.5078 | 0.1601 | 0.4857 |

The phage-prevalence baseline uses only each phage's infection frequency on
training strains, with Beta(1,1) smoothing for probability prediction.

This is a strong nuisance baseline: a future host-aware model must demonstrate
that it learns host-phage compatibility beyond merely identifying broad-host-
range phages.

## Offline broadest-phage union-coverage baseline

Phages are ranked using training strains only.

| k | Selected phages | Naive independence estimate | Held-out single-phage union coverage |
| ---: | --- | ---: | ---: |
| 1 | DIJ07_P2 | 67.74% | 61.54% |
| 2 | DIJ07_P2, 536_P7 | 88.03% | 73.08% |
| 3 | DIJ07_P2, 536_P7, 536_P9 | 95.46% | 73.08% |
| 4 | DIJ07_P2, 536_P7, 536_P9, DIJ07_P1 | 98.21% | 73.08% |
| 5 | DIJ07_P2, 536_P7, 536_P9, DIJ07_P1, LF73_P1 | 99.29% | 88.46% |

The independence formula becomes extremely optimistic for k=3 and k=4:
adding additional broad phages raises the mathematical independence estimate
while producing no additional test-strain coverage. This is evidence that
cocktail design must account for redundant host ranges and dependencies.

Coverage is calculated only when the test matrix has an explicit observed outcome
for every selected-phage/held-out-strain pair. Missing pairs are unknown, not
negative; evaluation fails closed rather than silently shrinking the denominator.

Union coverage is calculated from separate pair labels, not measured mixtures.
It cannot establish cocktail synergy, antagonism or efficacy.

These numbers are a fixed development benchmark, not a final estimate of model
performance. The source label's exact assay definition and mapping to
adsorption/productive-infection/resistance-emergence targets are not yet
verified. Later evaluation should include repeated group splits, confidence
intervals, and a separate unseen-phage protocol.
