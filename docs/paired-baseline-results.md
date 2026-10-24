# First additive host–phage baseline: results

**The prespecified cold-phage improvement criterion was not met.** The primary
pair model has AP 0.4761 versus 0.3985 for pooled host-only logistic, but the
paired phage-group bootstrap AP-gain interval includes zero:
**+0.0776, 95% interval [-0.1313, 0.2688]**. No sensitivity configuration establishes
a positive lower bound. These remain uncalibrated source-MLC research scores.

The model consumes features of both entities and can calculate scores for a phage
ID absent from fitting when its sequence-derived features are supplied. This
interface capability does not establish reliable unseen-phage performance.

## Protocol and input admission

The execution [protocol](../data/studies/paired_baseline.v2.json) binds the unchanged
v2 plan and frozen split manifests. All nine identity-threshold/seed schemes were
run with two prespecified observation-weighting policies, giving 18 scenarios.
Each scenario uses the same fit rows across its applicable models. No test metric
changed features, hyperparameters, thresholds, seeds or the selection rule.

All 276 input resources were re-admitted by the source lock. Missing MLC is
excluded, leaving 16,597 observations; the 65 conflicting exact sequence-pair
labels are retained. Raw genomes, labels, audit diagnostics and future outcomes
are not invented or filled. Published defense counts/surface categories are
inputs; the model does not run or interpret a defense annotation tool.

The [reuse review](../data/studies/paired_data_reuse_review.v1.json) records publisher
license declarations independently of software licensing. The original 146,786,594
byte Zenodo archive was downloaded and its published MD5 verified. Its 94 phage
genomes and four main analysis tables match the locked source files byte for byte.
All 177 selected host files match the Figshare version's file ID/size/MD5 metadata.
Both deposits declare CC BY 4.0. The separately unresolved GenoPHI-subset terms
are not declared cleared by its software license; fitted labels are the reconciled
primary MLC observations. Raw source files are not redistributed with these results.
The earlier record of an incomplete archive transfer remains a historical record.

## Representations and fitting

- Host: 137 defense-subtype counts plus O-type, H-type, LPS type and ABC capsule
  category, with explicit categorical missingness.
- Phage: 136 canonical DNA 4-mer frequencies. Contig boundaries are respected;
  ambiguous windows are skipped; no valid windows is unsupported, not an all-zero
  biological feature vector.
- Fit scaling/vocabulary: unique fit host IDs and unique fit phage assembly
  digests only. Evaluation entities do not affect preprocessing statistics.
- Pooled models: host-only, phage-only and concatenated additive logistic. C is
  selected from 0.01, 0.1 and 1 using cold-phage validation AP only. The primary
  scenario selected C=1 for all three; that is a selection result, not a new grid.
- Comparators: train-only global and per-phage prevalence, and the old per-phage
  host model refitted on the new fit rectangle. Known-phage comparators abstain
  on cold-phage/double-cold rows; frozen all-phage weights are not reused there.

The primary model was fit on 6,834 observations. The additive model represents
both entities but has no explicit host–phage interaction terms, RBP representation,
antidefense/prophage block or foundation-model embedding.

## Primary test results

Identity threshold 95%, bilateral coverage 85%, seed 42, equal observation weights:

| Regime | Rows | Pooled host AP | Phage 4-mer AP | Additive pair AP | Refit per-phage host AP |
| --- | ---: | ---: | ---: | ---: | ---: |
| Warm pair | 370 | 0.5469 | 0.4635 | 0.7754 | 0.4803 |
| Cold host | 2,210 | 0.3048 | 0.4901 | 0.4699 | 0.5468 |
| Cold phage | 1,633 | 0.3985 | 0.3381 | 0.4761 | Unsupported |
| Double cold | 476 | 0.3486 | 0.2937 | 0.3885 | Unsupported |

The primary cold-phage panel has ten held-out phage groups. Its global-prevalence
AP is 0.1923. Per-phage prevalence and the refitted known-phage host model have
zero support on that panel, rather than fabricated fallback scores.

| Regime | Pair-minus-pooled-host AP | Paired group bootstrap 95% interval |
| --- | ---: | --- |
| Warm pair | +0.2285 | [0.0578, 0.4183] |
| Cold host | +0.1651 | [0.1133, 0.2521] |
| Cold phage: primary | +0.0776 | [-0.1313, 0.2688] |
| Double cold | +0.0399 | [-0.1175, 0.2950] |

The secondary positive intervals compare against pooled host-only logistic, not
the strongest applicable comparator. On cold-host rows, per-phage prevalence AP
is 0.5066 and refitted per-phage host AP is 0.5468; both exceed the pair model's
0.4699. The new model therefore does not justify replacing the existing known-
phage baseline.

Ranking improvements also do not establish calibration. On double-cold rows,
the pair model has worse Brier/log loss (0.1899/0.6380) than pooled host-only
(0.1715/0.5601), despite its higher AP. Probability promotion was not performed.

## All prespecified sensitivities

Cold-phage AP gain over pooled host-only, equal observation weights:

| Identity | Seed | AP gain | Group bootstrap 95% interval |
| --- | ---: | ---: | --- |
| 95%: primary | 42 | +0.0776 | [-0.1313, 0.2688] |
| 95% | 43 | -0.0618 | [-0.1405, 0.0524] |
| 95% | 44 | -0.0537 | [-0.1199, 0.0188] |
| 90% | 42 | +0.1376 | [-0.1085, 0.2203] |
| 90% | 43 | -0.1151 | [-0.1511, -0.0124] |
| 90% | 44 | -0.0649 | [-0.1219, 0.0594] |
| 99% | 42 | +0.1563 | [-0.0902, 0.2386] |
| 99% | 43 | -0.0215 | [-0.1827, 0.1525] |
| 99% | 44 | +0.0227 | [-0.0750, 0.1192] |

Giving every exact sequence pair equal total weight leaves the conclusion
unchanged: no positive cold-phage lower bound in any of the nine schemes.
The primary weighted gain is +0.0859, interval [-0.1241, 0.2722]. All 18 contrasts
are retained in the saved summary; no favorable seed/weight policy is selected.
These sensitivity intervals are descriptive development comparisons, without a
multiplicity-adjusted claim from selecting one secondary result.

## Artifacts, verification and limits

The [result directory](../data/studies/results/paired_baseline.v2/) contains
aggregate metrics for every scenario, all C trials, all contrasts and a compressed
full run report. Primary development weights are saved separately in `models/`;
the original frozen host-baseline weights remain unchanged. Full local runs also save entity
features, models, row predictions, exclusions and provenance together.

An initial v1 numerical diagnostic found that batch BLAS scoring could split
identical feature rows by roughly 1e-16 and thereby alter rank ties. A failing
regression reproduced the issue. Scoring now accumulates features in a fixed
float64 order. All 18 unchanged planned scenarios were rerun under the v2 execution
record; the v1 plan, report and checkpoint are retained as **diagnostic, not validated
results**. No data, representation, parameter grid or selection rule was changed.
The primary conclusion remains unchanged.

Saved v2 weights/encoders reproduce the stored primary scores bit for bit, including
across batch shapes. Exact host-only/phage-only ties were checked in every regime
and partition of all 18 scenarios. Run protocol, source locks, code versions, outputs and
split hashes were verified. Tests cover strand/contig handling, held-out statistics,
sequence-pair weights, unsupported outputs, test-label independence and bootstrap
ties/crossed axes. Nonconvergence is a failure, with no silent alternate model.

The scores predict source-coded MLC positivity. They do not establish productive
infection, cocktail effects, clinical benefit, calibrated probabilities, OOD
rejection or independent release validation. The current file-based exchange
still has no serving adapter for these research weights. Model selection and
scientific promotion remain separate from exchange integrity.

To reproduce in the pinned numerical environment:

```bash
pip install -e ".[dev]" -c data/studies/host_baseline_environment.v1.txt
python scripts/run_paired_baseline.py --fetch --outdir artifacts/paired_baseline_v1
```

Use a new output directory. The fixed manifests and execution protocol must match
their hashes. This comparison does not test whether explicit interaction models,
RBP/protein features or pretrained representations would perform better.
