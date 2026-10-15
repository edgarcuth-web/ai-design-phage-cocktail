# Data provenance

Large raw datasets are intentionally not committed to this repository.

## Step 1 dataset: GenoPHI E. coli subset

Source project: GenoPHI by Noonan et al.

Pinned upstream commit:

    f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f

Upstream file:

    data/interaction_matrices/ecoli_interaction_matrix_subset.csv

Download it reproducibly with:

    python scripts/fetch_genophi.py

The upstream schema is:

    strain,phage,interaction

The project loader normalizes it to:

    bacterium,phage,infection

Expected statistics for this pinned file:

- 16,638 measured interactions
- 177 bacterial strains
- 94 phages
- 3,813 positive interactions
- 12,825 negative interactions
- 22.9% positive rate
- 100% of the 177 x 94 strain-phage matrix observed
- no duplicate strain-phage pairs

Publication:

Noonan AJC et al. "Phylogeny-agnostic strain-level prediction of phage-host
interactions from genomes using machine learning." Nature Microbiology (2026).

The processed publication data are also archived on Zenodo:
DOI 10.5281/zenodo.21696332.
