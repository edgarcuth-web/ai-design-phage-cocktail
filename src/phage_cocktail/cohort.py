"""Observed-label curation and genomic grouping for a pinned research cohort."""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

from .data import validate_interactions


def host_aliases(names) -> dict[str, str]:
    """Only the reviewed case/hyphen transformation is allowed; no fuzzy joins."""
    result = {}
    for name in names:
        name = str(name)
        key = name.replace("-", "").upper()
        if not key or key in result:
            raise ValueError(f"ambiguous host alias: {name}")
        result[key] = name
    return result


def curate_mlc(pairs: pd.DataFrame, published: pd.DataFrame):
    """Check all candidate labels; return known rows, explicit exclusions and IDs."""
    validate_interactions(pairs)
    if published.columns.duplicated().any():
        raise ValueError("duplicate published phage columns")
    lookup = host_aliases(published.index)
    aliases = {}
    for host in sorted(pairs.bacterium.astype(str).unique()):
        key = host.replace("-", "").upper()
        if key not in lookup:
            raise ValueError(f"host missing from primary MLC source: {host}")
        aliases[host] = lookup[key]
    if len(set(aliases.values())) != len(aliases):
        raise ValueError("ambiguous cohort host aliases")
    observed, excluded = [], []
    for row in pairs.itertuples(index=False):
        if row.phage not in published:
            raise ValueError(f"phage missing from primary MLC source: {row.phage}")
        value = published.loc[aliases[row.bacterium], row.phage]
        common = {"bacterium": row.bacterium, "source_bacterium": aliases[row.bacterium],
                  "phage": row.phage, "source_binary": int(row.infection)}
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
    return (pd.DataFrame(observed, columns=columns),
            pd.DataFrame(excluded, columns=[*columns[:4], "reason"]), aliases)


def canonical_assembly_digest(path: str | Path) -> str:
    """Hash contig content independent of FASTA headers, order and strand."""
    complements = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    allowed = set("ACGTRYSWKMBDHVN")
    contigs, lines, seen_header = [], [], False

    def finish():
        sequence = "".join(lines).upper()
        if not sequence or set(sequence) - allowed:
            raise ValueError(f"invalid or empty DNA contig in {path}")
        canonical = min(sequence, sequence.translate(complements)[::-1])
        contigs.append(sha256(canonical.encode("ascii")).hexdigest())

    with Path(path).open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    finish()
                seen_header, lines = True, []
            else:
                if not seen_header:
                    raise ValueError(f"DNA FASTA requires a header: {path}")
                lines.append(line)
    if not seen_header:
        raise ValueError(f"empty DNA FASTA: {path}")
    finish()
    return sha256("\n".join(sorted(contigs)).encode("ascii")).hexdigest()


def group_hosts(hosts, distances: pd.DataFrame, assembly_digests: dict, *, threshold: float):
    """Union exact assembly matches and published near-neighbours transitively."""
    hosts = sorted(set(map(str, hosts)))
    if threshold <= 0 or not np.isfinite(threshold):
        raise ValueError("core distance threshold must be positive and finite")
    if set(hosts) - set(assembly_digests):
        raise ValueError("every host requires a verified assembly digest")
    row_aliases, col_aliases = host_aliases(distances.index), host_aliases(distances.columns)
    try:
        rows = [row_aliases[h.replace("-", "").upper()] for h in hosts]
        cols = [col_aliases[h.replace("-", "").upper()] for h in hosts]
    except KeyError as exc:
        raise ValueError(f"missing core distance host: {exc}") from exc
    values = distances.loc[rows, cols].to_numpy(dtype=float)
    if (not np.isfinite(values).all() or (values < 0).any()
            or not np.allclose(values, values.T, atol=1e-10, rtol=0)
            or not np.allclose(np.diag(values), 0, atol=1e-10, rtol=0)):
        raise ValueError("core distances must be finite, non-negative, symmetric with zero diagonal")
    parent = {h: h for h in hosts}

    def find(h):
        while parent[h] != h:
            parent[h] = parent[parent[h]]
            h = parent[h]
        return h

    for i, left in enumerate(hosts):
        for j, right in enumerate(hosts[:i]):
            if values[i, j] < threshold or assembly_digests[left] == assembly_digests[right]:
                a, b = find(left), find(right)
                parent[max(a, b)] = min(a, b)
    return {h: find(h) for h in hosts}
