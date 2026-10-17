# Data provenance

Raw and processed datasets are intentionally not committed. The pinned source,
checksum, label semantics and expected summary for the first dataset are recorded
in the machine-readable registry at
[`registry/genophi_ecoli_subset.v1.json`](registry/genophi_ecoli_subset.v1.json).

## Step 1 dataset: GenoPHI E. coli subset

Source project: GenoPHI by Noonan et al.

Pinned upstream commit:

    f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f

Upstream file:

    data/interaction_matrices/ecoli_interaction_matrix_subset.csv

Download it reproducibly with:

    python scripts/fetch_genophi.py

The download script verifies SHA-256 before writing the file. Analysis scripts
verify the checksum and expected summary again before using the dataset.

The upstream schema is:

    strain,phage,interaction

The project loader normalizes it to:

    bacterium,phage,infection

This is a column-name normalization, not a claim that the source label directly
measures adsorption, productive infection, or resistance emergence. The exact
source assay definition and mapping to those mechanistic targets remain
unverified and must be reviewed before the label is used for those tasks.

Expected statistics for this pinned file:

- 16,638 observed interactions
- 177 bacterial strains
- 94 phages
- 3,813 source-positive labels
- 12,825 source-negative labels
- 22.9% positive rate
- 100% of the 177 x 94 strain-phage matrix observed
- no duplicate strain-phage pairs

**Missing-pair rule:** a pair absent from a dataset is unknown/unobserved, not
negative. Do not fill missing pairs with zero. The source subset happens to be a
complete matrix; future datasets may not be.

License status has not yet been independently verified. Do not redistribute the
downloaded data until the upstream data terms have been checked.

Publication:

Noonan AJC et al. "Phylogeny-agnostic strain-level prediction of phage-host
interactions from genomes using machine learning." Nature Microbiology (2026).

The processed publication data are also archived on Zenodo:
DOI 10.5281/zenodo.21696332.
