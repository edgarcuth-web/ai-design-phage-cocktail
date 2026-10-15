# Step 2 benchmark

This benchmark uses the pinned GenoPHI E. coli interaction matrix and a
deterministic bacterial-strain split with seed 42.

## Split

| Partition | Strains | Interactions | Positive rate |
| --- | ---: | ---: | ---: |
| Train | 124 | 11,656 | 23.14% |
| Validation | 27 | 2,538 | 18.72% |
| Test | 26 | 2,444 | 26.23% |

No bacterial strain occurs in more than one partition.

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

## Broadest-phage cocktail baseline

Phages are ranked using training strains only.

| k | Selected phages | Naive independence estimate | Actual held-out coverage |
| ---: | --- | ---: | ---: |
| 1 | DIJ07_P2 | 67.74% | 61.54% |
| 2 | DIJ07_P2, 536_P7 | 88.03% | 73.08% |
| 3 | DIJ07_P2, 536_P7, 536_P9 | 95.46% | 73.08% |
| 4 | DIJ07_P2, 536_P7, 536_P9, DIJ07_P1 | 98.21% | 73.08% |
| 5 | DIJ07_P2, 536_P7, 536_P9, DIJ07_P1, LF73_P1 | 99.29% | 88.46% |

The independence formula becomes extremely optimistic for k=3 and k=4:
adding additional broad phages raises the mathematical independence estimate
while producing no additional test-strain coverage. This is direct evidence
that cocktail design must account for redundant host ranges and dependencies.

These numbers are a fixed development benchmark, not a final estimate of model
performance. Later evaluation should include repeated group splits or grouped
cross-validation and confidence intervals.
