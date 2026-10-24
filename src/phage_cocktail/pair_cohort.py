"""Audit paired input identities without fitting models or inferring labels."""
from __future__ import annotations


from .data import validate_interactions


def audit_pair_identities(observed, host_digests, phage_digests):
    """Group exact phage assemblies; retain conflicting sequence-pair labels.

    Digests must come from admitted sequences. Matching sequence content alone
    does not establish matching physical material, assay context or phenotype.
    """
    validate_interactions(observed)
    for subject, digests in (("bacterium", host_digests), ("phage", phage_digests)):
        if set(observed[subject]) != set(digests) or any(not d for d in digests.values()):
            raise ValueError(f"{subject} genome IDs must match the observed cohort exactly")
    representatives = {}
    for phage in sorted(phage_digests):
        representatives.setdefault(phage_digests[phage], phage)
    groups = {phage: representatives[digest] for phage, digest in phage_digests.items()}
    rows = observed.copy()
    rows["host_assembly_sha256"] = rows.bacterium.map(host_digests)
    rows["phage_assembly_sha256"] = rows.phage.map(phage_digests)
    rows["exact_phage_group"] = rows.phage.map(groups)
    keys = ["host_assembly_sha256", "phage_assembly_sha256"]
    conflicting = rows.groupby(keys).infection.transform("nunique") > 1
    conflicts = rows.loc[conflicting].sort_values(["bacterium", "phage"]).reset_index(drop=True)
    return groups, conflicts
