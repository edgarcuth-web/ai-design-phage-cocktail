# Endpoint semantics

This note defines interpretation and review requirements for the current exchange.
It adds no wire fields, endpoint mappings or supported model targets. The same note
is maintained in SDL and `phage-host-interaction-predictor`.

## Readouts and biological questions

An endpoint names the question and decision rule; an assay supplies observations.
A result is meaningful only with its source definition and assay context. The
following distinctions guide review; they are not claims that these measurements
are already collected or predicted by either repository.

| Endpoint / readout | Interpretation to preserve | Do not infer automatically |
| --- | --- | --- |
| Source-coded binary label, including the audited `MLC > 0` target | Positive under the pinned source's coding/curation rule | Adsorption, productive infection, or general susceptibility |
| Spotting / clearing | Clearing under the stated assay conditions and scoring rule | Plaque formation or production of infectious progeny |
| Plaque formation | Plaques under the stated assay conditions and scoring rule | An interchangeable label for spotting, liquid growth suppression, or every productive-infection assay |
| Adsorption | The explicitly measured binding/adsorption endpoint | Genome entry or successful intracellular replication |
| Genome entry | The explicitly measured entry endpoint | Completion of replication or infectious progeny production |
| Productive infection | An explicitly defined endpoint supported by evidence of replication and infectious progeny production | A direct synonym for any clearing or growth-inhibition observation |
| Growth suppression | The declared change in bacterial growth over a specified observation period and comparator | Infectious progeny production or durable suppression of resistance |
| Resistance emergence | A defined change in susceptibility over a declared time horizon, exposure history and assessment | A label derived from a single initial positive/negative interaction |

The same experiment may support more than one interpretation after explicit
assessment. The software must retain the observations and the reasoning for each
label rather than assuming equivalence. These are review distinctions, not a
prescribed wet-lab procedure.

## Current labels and missingness

The frozen historical endpoint is
`genophi-ecoli-subset:source-interaction:v1`, source column `interaction`.
The canonical table's `infection` name is a rename of a source-coded label,
not a mechanistic claim. The later host-baseline study reconciles observed source
readings with `MLC > 0` and excludes 41 source-coded zeros with missing primary
readings. This changes the evaluated cohort; it does not retroactively verify or
replace the frozen historical endpoint.

The predictor's `docs/first-model-task.md`, `docs/host-baseline-results.md` and
`data/studies/genophi_host_baseline.v1.json` record the audit and its remaining
limits. The frozen model's `uncalibrated_source_mlc_positive` scores refer to that
audited target. Row-specific assay context and mechanistic mappings remain
unverified. A different endpoint/curation rule needs an explicitly reviewed
version and provenance, not an edit to an already frozen descriptor.

Missing observations, absent annotations, unresolved predictions, invalid QC and
untested pairs are not negative labels. Model abstention describes a model's
output state; it does not describe an observed biological outcome.

## Existing software boundaries

| Boundary | Current behavior | Scientific review still required |
| --- | --- | --- |
| ModelInputBundle v1 | Binds a source endpoint, exact pair, assay context and inputs by content hash | Source protocol, decision rule, missingness and relevance to the requested task |
| PairPrediction v1 | Binds the exact input bundle; distinguishes unresolved, score and probability | Model target matches the requested source endpoint; applicability and evaluation evidence |
| Source probability admission | Requires a declared verified assay definition, calibration digest and containing uncertainty interval | Whether those claims and artifacts are valid; calibration alone does not map the endpoint to a mechanism |
| Legacy ModelResult / InteractionPrediction | Represents adsorption, productive infection and resistance emergence separately | Target-specific label evidence, calibration and assay comparability |
| CandidateSet screening policy | Accepts InteractionPrediction; ranks productive-infection lower bounds and retains an unresolved exploration queue | Whether the supplied mechanistic estimates are scientifically suitable for that context |
| LearningObservation | Requires an explicit assessor target/label and exact physical phenotype lineage | Whether the observation actually supports the declared target |

A structurally admissible source probability is still a probability for its source
endpoint. Neither renaming fields nor supplying calibration provenance admits it
to the productive-infection policy. There is currently no reviewed mapping or
source-endpoint candidate-selection policy. Current host-baseline weights are not
served through this bundle interface and are not admitted to SDL decisions.

Changing an endpoint or assay context changes the bundle hash and requires a new
export/result. Hash agreement establishes binding, not a valid scientific mapping.
A model implementation must check target compatibility before inference; the
generic exchange validator does not inspect the training labels behind weights.

## Requirements before admitting a new target or mapping

Record the source observation, positive/negative decision rule, exclusions and
unknowns, exact host/phage identities, applicable conditions and time horizon,
and supporting evidence/review. Preserve original readouts alongside any derived
label and document the mapping's limits. The existing endpoint definition,
evidence/limitations and assay-context fields can carry reviewed descriptions;
they do not provide dedicated structured fields for every item in this checklist.

Then validate the target on an appropriate evaluation cohort, review calibration
and uncertainty if probabilities are claimed, and explicitly define the receiving
decision policy. Do not infer productive-infection labels from growth or clearing,
resistance labels from initial susceptibility, or individual-member labels from a
cocktail outcome. Recording new evidence does not automatically retrain or release
a model.
