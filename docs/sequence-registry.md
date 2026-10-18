# Local sequence registry v1

## Purpose and boundary

The resource resolver validates exact HostState and PhageProfile JSON records.
This layer resolves the raw FASTA assets identified by those records and verifies
both the file bytes and the normalized sequence content.

Sequence resolution is an integrity and identity check only. It does not confirm
that the source database annotation is correct, that a sequence belongs to the
claimed organism, that the isolate is pure, or that any biological feature or
interaction prediction is valid. It does not fetch remote content, derive
features, train models, or run inference.

## Two different SHA-256 digests

Each registered FASTA has two distinct fingerprints:

- `file_sha256`: SHA-256 of the exact original file bytes, including headers,
  line endings, wrapping, case, and any UTF-8 byte-order mark.
- `sequence_sha256`: SHA-256 of normalized contigs, matching the existing lab
  `phage_lab/biology.py::read_fasta` convention.

The normalized digest is computed as follows:

1. Decode UTF-8, accepting an optional UTF-8 BOM, and use universal newline
   handling.
2. Ignore blank lines and FASTA header text. A header is a non-empty line
   beginning with `>` after surrounding whitespace is stripped.
3. For each header, append one NUL byte (`0x00`) to the digest input.
4. For each sequence line, strip surrounding whitespace, convert bases to
   uppercase, validate against `ACGTRYSWKMBDHVN`, then append the ASCII bases.
5. The normalized base count is the sum of bases across all contigs.

In other words, for contigs `ACGT` and `NRY`, the hashed bytes are
`NUL + ACGT + NUL + NRY`. Headers, whitespace, line wrapping and letter case do
not affect `sequence_sha256`, but they do affect `file_sha256`. Each contig
boundary is retained, so different contig segmentation does affect the digest.

The current limits match the lab's genome reader: 20 million total bases and
40,000,000 original file bytes. Empty contigs, empty headers, sequence before a
header, invalid symbols, malformed UTF-8, and empty FASTA files are rejected.
FASTA headers are not used as an identity source.

## Manifest format

A sequence registry is a JSON object with exactly these top-level fields:

- `schema_version`: integer `1`;
- `registry_id`: a non-empty local identifier;
- `entries`: a list of local FASTA mappings.

The machine-readable contract is
[`sequence-registry-v1.schema.json`](../schemas/sequence-registry-v1.schema.json).

Every entry requires:

- `kind`: `host_genome` or `phage_sequence`;
- `sequence_sha256`: expected normalized sequence digest;
- `file_sha256`: expected original file-byte digest;
- `path`: POSIX-style path relative to the sequence registry file;
- `status`: `available` or `withdrawn`.

A phage-sequence entry may also include `accession` and
`source_snapshot_sha256`. When the PhageProfile names an accession, the
registered accession must match. If the profile also names a retrieved-source
snapshot SHA-256, the entry must carry the same snapshot digest. Host-genome
entries cannot declare phage accession or snapshot fields.

Host-genome digests must be unique; phage accession plus source-snapshot digest
must identify at most one entry; file paths must be unique. Absolute paths,
`.`/`..` components, backslashes, missing files, and symlinks that escape the
registry directory are rejected. The resolver only uses paths from the locally
managed registry.

## Resolution rules

The resolver first validates the ModelRequest and resolves its exact HostState
and PhageProfile records using the resource registry. For retrieved profiles,
the sequence accession must agree with provenance, and snapshot version/hash must
be supplied together. It then:

1. resolves the HostState's `host_genome_sha256` as a `host_genome` asset;
2. requires the PhageProfile to have a sequence reference with an explicit
   `sha256`;
3. resolves that fingerprint as a `phage_sequence` asset;
4. checks an optional profile accession against the registered accession;
5. verifies original file bytes and recomputes the normalized sequence digest
   from the same bounded byte buffer;
6. checks the profile's optional sequence length against the observed base count.

A phage reference with a normalized sequence SHA-256 resolves by content digest,
and any accession/snapshot fields present must match. If the profile has only an
accession, resolution is allowed only for a `retrieved` profile with a pinned
`provenance.snapshot_sha256`; the registry must bind that exact accession and
snapshot digest to one FASTA asset. The resolver checks the profile's accession
against its retrieved provenance. An accession alone is not sufficient.

The resolver does not independently re-hash the full upstream snapshot. The
trusted ingestion process that creates the sequence registry must verify the
snapshot against its recorded digest before extracting the accession-specific
FASTA. The registry records that provenance binding alongside the extracted
file's own byte and normalized sequence hashes. A missing PhageProfile sequence
reference is not treated as a negative or a successful resolution.

The audit summary reports both digests, normalized base counts, contig counts,
accession when present, and registry-relative paths. It never includes the raw
sequence content. It explicitly records that features have not been derived and
biological validity has not been established.

## Usage

After installing the project, prepare a sequence registry only from reviewed
FASTA assets. The values for `sequence_sha256` must use the normalized digest
convention above; `file_sha256` must hash the exact file bytes.

Run:

```bash
phage-model-interface resolve-request-sequences \
  path/to/model-request.json \
  --registry path/to/resource-registry.v1.json \
  --sequence-registry path/to/sequence-registry.v1.json
```

The command prints an audit summary for the request, both record registries and
the two resolved sequence assets. It does not copy or print FASTA sequences.

## Next boundary

The provenance manifest is now defined in [FeatureArtifact v1](feature-artifact.md). It records normalized and raw input digests, tool/version/parameters, reference snapshots, output digests and limitations. No analyzer is wired into this contract yet. The next step is a narrow adapter around one reviewed real tool, after checking the intended feature class, pinned database release and license.
