# Research weights

`host_logistic_development.v1.json` is the frozen first host-feature baseline,
trained on the audited source MLC cohort. It contains a training-only encoder,
one model per supported phage, explicit single-class fallback and exact input,
code and environment provenance. `phage_cocktail.host_baseline.predict_host_model`
loads these JSON weights without executing serialized Python objects.

Inputs are the published defense counts and four surface/capsule annotations
listed in the encoder, plus known phage identity. A raw FASTA or an ID/hash alone
does not supply these required feature values. Unknown phages are unsupported;
missing defense counts fail rather than becoming observed absence.

The model did not meet the predeclared criterion for beating phage prevalence.
Its outputs remain uncalibrated source-label scores. It is not admitted into
SDL candidate decisions. See [the full study report](../docs/host-baseline-results.md)
for evaluation, limitations, hashes and reproduction.

## Paired development weights

`paired_development.v2.json.gz` contains the primary 95%/seed-42/equal-observation
study models: one train-only entity encoder, selected host-only/phage-only/additive
logistic weights, prevalence comparators and the refitted known-phage host model.
It is JSON compressed with gzip, with no serialized Python code. The old host
weights above are a different frozen study and are preserved.

The host input requires the named published defense counts and four surface
categories. The phage input requires 136 canonical 4-mer frequencies extracted
with `canonical_fourmers`; a raw ID/hash alone does not supply those values.
`PairEncoder.from_dict` and `predict_logistic` load the three pooled models from
the saved `encoder` and `models` dictionaries. A new phage ID can be scored with
supplied features, but cold-phage improvement was not demonstrated. This is an
additive compositional baseline, with no explicit pair interaction terms.

The primary cold-phage pair AP is 0.4761 versus 0.3985 for pooled host-only,
AP-gain interval [-0.1313, 0.2688]. All 18 prespecified threshold/seed/weight
configurations fail the positive-lower-bound criterion. Weights are not calibrated,
promoted, served through PairPrediction, or admitted to SDL decisions. See the
[results](../docs/paired-baseline-results.md) and
[execution protocol](../data/studies/paired_baseline.v2.json).

`paired_development.v1.json.gz` is retained only for the documented batch-score
tie diagnostic. Its associated results are not the validated comparison; use v2.
