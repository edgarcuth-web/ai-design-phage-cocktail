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
