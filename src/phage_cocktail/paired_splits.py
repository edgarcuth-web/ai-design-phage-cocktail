"""Operational genome-neighbour grouping and crossed development partitions."""
from collections import defaultdict
from hashlib import sha256

import pandas as pd

from .cohort import DisjointSet
from .data import validate_interactions
from .split import assign_partitions


def neighbour_groups(digests, blocks, *, identity, coverage):
    """Union exact assemblies and bilateral non-overlapping alignment coverage.

    These components control a specified similarity boundary, not taxonomy or
    functional relatedness. Input alignments must be base-level, not sketches.
    """
    if not 0 < identity <= 1 or not 0 < coverage <= 1:
        raise ValueError("identity and coverage must lie in (0,1]")
    components, edges = DisjointSet(digests), set()

    def union(a, b):
        if a != b:
            edges.add(tuple(sorted((a, b))))
            components.union(a, b)

    exact = {}
    for name in sorted(digests):
        if not isinstance(digests[name], str) or not digests[name]:
            raise ValueError("every genome needs an admitted assembly digest")
        digest = digests[name]
        if digest in exact:
            union(name, exact[digest])
        exact[digest] = name
    grouped = defaultdict(list)
    for row in blocks:
        q, t = row["query"], row["target"]
        if q not in digests or t not in digests:
            raise ValueError("alignment has an unknown genome")
        if not (
            0 <= row["qs"] < row["qe"] <= row["qlen"]
            and 0 <= row["ts"] < row["te"] <= row["tlen"]
            and 0 <= row["matches"] <= row["block"]
            and row["block"] > 0
        ):
            raise ValueError("invalid alignment coordinates or identity counts")
        if q != t and row["matches"] / row["block"] >= identity:
            grouped[q, t].append(row)
    for (q, t), hits in sorted(grouped.items()):
        if len({(h["qlen"], h["tlen"]) for h in hits}) != 1:
            raise ValueError("inconsistent genome lengths in alignments")
        # Greedily keep the strongest blocks that overlap no kept block on either genome.
        retained = []
        best_first = sorted(
            hits,
            key=lambda h: (-h["matches"], -h["block"], h["qs"], h["qe"], h["ts"], h["te"], h["strand"]),
        )
        for row in best_first:
            overlaps = any(
                (row["qs"] < old["qe"] and old["qs"] < row["qe"])
                or (row["ts"] < old["te"] and old["ts"] < row["te"])
                for old in retained
            )
            if not overlaps:
                retained.append(row)
        qcov = sum(h["qe"] - h["qs"] for h in retained) / hits[0]["qlen"]
        tcov = sum(h["te"] - h["ts"] for h in retained) / hits[0]["tlen"]
        if min(qcov, tcov) >= coverage:
            union(q, t)
    return components.labels(), sorted(edges)


def parse_paf(path, lengths):
    """Require base-level CIGAR and lengths that match the admitted assemblies."""
    blocks = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 13 or not any(tag.startswith("cg:Z:") for tag in fields[12:]):
            raise ValueError("PAF requires base-level CIGAR, not approximate mapping")
        q, t = fields[0], fields[5]
        if q not in lengths or t not in lengths:
            raise ValueError("PAF has an unknown genome")
        qlen, qs, qe = map(int, fields[1:4])
        tlen, ts, te, matches, block = map(int, fields[6:11])
        if qlen != lengths[q] or tlen != lengths[t] or fields[4] not in {"+", "-"}:
            raise ValueError("PAF genome lengths/strand differ from inputs")
        blocks.append(
            dict(
                query=q, target=t, qlen=qlen, qs=qs, qe=qe, tlen=tlen, ts=ts, te=te,
                matches=matches, block=block, strand=fields[4],
            )
        )
    return blocks


def _partitions(groups, seed):
    """Partition genomes by assigning their whole groups."""
    if any(not isinstance(group, str) or not group for group in groups.values()):
        raise ValueError("every genome requires a non-empty group")
    by_group = assign_partitions(groups.values(), seed=seed)
    if set(by_group.values()) != {"train", "validation", "test"}:
        raise ValueError("insufficient groups for three non-empty partitions")
    return {name: by_group[group] for name, group in groups.items()}


def assign_paired_rows(observed, host_groups, phage_groups, host_digests, phage_digests, *, seed):
    """Allocate every observed row once; withhold both crossed test axes from fit."""
    validate_interactions(observed)
    for subject, groups, digests in (
        ("host", host_groups, host_digests),
        ("phage", phage_groups, phage_digests),
    ):
        column = "bacterium" if subject == "host" else "phage"
        if set(observed[column]) != set(groups) or set(groups) != set(digests):
            raise ValueError(f"{subject} grouping/digests must cover the exact cohort")
        exact = {}
        for name, digest in digests.items():
            if not isinstance(digest, str) or not digest:
                raise ValueError("every genome requires a non-empty digest")
            if digest in exact and groups[name] != exact[digest]:
                raise ValueError(f"identical {subject} genomes cannot cross groups")
            exact[digest] = groups[name]
    hosts, phages = _partitions(host_groups, seed), _partitions(phage_groups, seed)
    rows = []
    for row in observed.itertuples(index=False):
        h, p = row.bacterium, row.phage
        hp, pp = hosts[h], phages[p]
        if hp == pp == "train":
            # Hold out ~5% + 5% of seen-host/seen-phage pairs by a stable hash of the sequence pair.
            key = f"{seed}:{host_digests[h]}:{phage_digests[p]}".encode()
            bucket = int(sha256(key).hexdigest(), 16) % 100
            partition = "validation" if bucket < 5 else "test" if bucket < 10 else "fit"
            regime = "fit" if partition == "fit" else "warm_pair"
        elif hp == "train":
            regime, partition = "cold_phage", pp
        elif pp == "train":
            regime, partition = "cold_host", hp
        elif hp == pp:
            regime, partition = "double_cold", hp
        else:
            regime, partition = "unused", "unused"
        rows.append(
            dict(
                bacterium=h,
                phage=p,
                host_group=host_groups[h],
                phage_group=phage_groups[p],
                host_partition=hp,
                phage_partition=pp,
                regime=regime,
                partition=partition,
            )
        )
    assigned = pd.DataFrame(rows).sort_values(["bacterium", "phage"]).reset_index(drop=True)
    fit = assigned[assigned.partition == "fit"]
    for column, partitions in (("bacterium", hosts), ("phage", phages)):
        intended = {name for name, part in partitions.items() if part == "train"}
        if set(fit[column]) != intended:
            raise ValueError(f"intended seen {column} disappeared from fit")
    fit_labels = observed.merge(
        fit[["bacterium", "phage"]], on=["bacterium", "phage"], validate="one_to_one"
    )
    if fit_labels.infection.nunique() != 2:
        raise ValueError("fit rectangle must have both observed classes")
    return assigned
