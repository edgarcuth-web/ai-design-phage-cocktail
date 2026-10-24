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
measures adsorption, productive infection, or resistance emergence. The
[later MLC audit](../docs/host-baseline-results.md) establishes exact binary
agreement on known primary readings and identifies 41 source zeros with missing
MLC. Mapping the target to mechanistic outcomes remains unverified.

Expected statistics for this pinned file:

- 16,638 source-coded rows
- 177 bacterial strains
- 94 phages
- 3,813 source-positive labels
- 12,825 source-negative labels
- 22.9% positive rate
- 100% of the 177 x 94 strain-phage matrix syntactically populated
- no duplicate strain-phage pairs

**Missing-pair rule:** a pair absent from a dataset is unknown/unobserved, not
negative. Do not fill missing pairs with zero. The source CSV is syntactically
complete, but 41 primary readings are missing. The new study admits 16,597
observed rows with 12,784 explicit negatives.

License status has not yet been independently verified. Do not redistribute the
downloaded data until the upstream data terms have been checked.

The later [paired-study review](studies/paired_data_reuse_review.v1.json) separately
records the primary Zenodo/Figshare CC BY 4.0 declarations and verifies archive
byte correspondence for 94 phage genomes and four analysis tables, plus selected
host-file metadata. This does not change the frozen GenoPHI registry's status or
treat a software license as a separate dataset license. Paired-study targets are
reconciled primary MLC observations.

Publication:

Noonan AJC et al. "Phylogeny-agnostic strain-level prediction of phage-host
interactions from genomes using machine learning." Nature Microbiology (2026).

The processed publication data are also archived on Zenodo:
DOI 10.5281/zenodo.21696332.
