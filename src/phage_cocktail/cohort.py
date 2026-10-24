"""Observed-label curation and genomic grouping for a pinned research cohort."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

from .data import validate_interactions

DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def host_key(name) -> str:
    """The single reviewed host-name normalization: drop hyphens, ignore case."""
    return str(name).replace("-", "").upper()


def host_aliases(names) -> dict[str, str]:
    """Map normalized keys to source names; reject collisions instead of fuzzy joins."""
    result = {}
    for name in names:
        name = str(name)
        key = host_key(name)
        if not key or key in result:
            raise ValueError(f"ambiguous host alias: {name}")
        result[key] = name
    return result


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def parse_dna_fasta(text: str, source="FASTA") -> list[str]:
    """Return upper-case contigs; every contig must be non-empty IUPAC DNA."""
    contigs, pieces = [], None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if pieces is not None:
                contigs.append("".join(pieces).upper())
            pieces = []
        elif pieces is None:
            raise ValueError(f"DNA FASTA requires a header: {source}")
        else:
            pieces.append(line)
    if pieces is None:
        raise ValueError(f"empty DNA FASTA: {source}")
    contigs.append("".join(pieces).upper())
    if any(not contig or set(contig) - DNA_ALPHABET for contig in contigs):
        raise ValueError(f"invalid or empty DNA contig in {source}")
    return contigs


def canonical_assembly_digest(path: str | Path) -> str:
    """Hash contig content independent of FASTA headers, order and strand."""
    contigs = parse_dna_fasta(Path(path).read_text(encoding="ascii"), path)
    digests = [
        sha256(min(contig, reverse_complement(contig)).encode("ascii")).hexdigest()
        for contig in contigs
    ]
    return sha256("\n".join(sorted(digests)).encode("ascii")).hexdigest()


class DisjointSet:
    """Union-find whose component label is the lexicographically smallest member."""

    def __init__(self, names):
        self.parent = {name: name for name in names}

    def find(self, name):
        parent = self.parent
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    def union(self, a, b):
        left, right = self.find(a), self.find(b)
        self.parent[max(left, right)] = min(left, right)

    def labels(self) -> dict:
        return {name: self.find(name) for name in sorted(self.parent)}


def curate_mlc(pairs: pd.DataFrame, published: pd.DataFrame):
    """Check all candidate labels; return known rows, explicit exclusions and IDs."""
    validate_interactions(pairs)
    if published.columns.duplicated().any():
        raise ValueError("duplicate published phage columns")
    lookup = host_aliases(published.index)
    aliases = {}
    for host in sorted(pairs.bacterium.astype(str).unique()):
        if host_key(host) not in lookup:
            raise ValueError(f"host missing from primary MLC source: {host}")
        aliases[host] = lookup[host_key(host)]
    if len(set(aliases.values())) != len(aliases):
        raise ValueError("ambiguous cohort host aliases")

    observed, excluded = [], []
    for row in pairs.itertuples(index=False):
        if row.phage not in published:
            raise ValueError(f"phage missing from primary MLC source: {row.phage}")
        value = published.loc[aliases[row.bacterium], row.phage]
        common = {
            "bacterium": row.bacterium,
            "source_bacterium": aliases[row.bacterium],
            "phage": row.phage,
            "source_binary": int(row.infection),
        }
        if pd.isna(value):
            excluded.append({**common, "reason": "source_mlc_missing"})
            continue
        if value not in {0, 1, 2, 3, 4}:
            raise ValueError(f"invalid MLC score: {value}")
        label = int(value > 0)
        if label != row.infection:
            raise ValueError(f"source labels disagree: {row.bacterium}/{row.phage}")
        observed.append({**common, "source_mlc": int(value), "infection": label})

    columns = ["bacterium", "source_bacterium", "phage", "source_binary", "source_mlc", "infection"]
    return (
        pd.DataFrame(observed, columns=columns),
        pd.DataFrame(excluded, columns=[*columns[:4], "reason"]),
        aliases,
    )


def group_hosts(hosts, distances: pd.DataFrame, assembly_digests: dict, *, threshold: float):
    """Union exact assembly matches and published near-neighbours transitively."""
    hosts = sorted(set(map(str, hosts)))
    if threshold <= 0 or not np.isfinite(threshold):
        raise ValueError("core distance threshold must be positive and finite")
    if set(hosts) - set(assembly_digests):
        raise ValueError("every host requires a verified assembly digest")
    row_aliases, col_aliases = host_aliases(distances.index), host_aliases(distances.columns)
    try:
        rows = [row_aliases[host_key(h)] for h in hosts]
        cols = [col_aliases[host_key(h)] for h in hosts]
    except KeyError as exc:
        raise ValueError(f"missing core distance host: {exc}") from exc
    values = distances.loc[rows, cols].to_numpy(dtype=float)
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.allclose(values, values.T, atol=1e-10, rtol=0)
        or not np.allclose(np.diag(values), 0, atol=1e-10, rtol=0)
    ):
        raise ValueError("core distances must be finite, non-negative, symmetric with zero diagonal")

    components = DisjointSet(hosts)
    for i, left in enumerate(hosts):
        for j, right in enumerate(hosts[:i]):
            if values[i, j] < threshold or assembly_digests[left] == assembly_digests[right]:
                components.union(left, right)
    return components.labels()
