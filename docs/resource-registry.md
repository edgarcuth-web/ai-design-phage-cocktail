# Local resource registry v1

## Purpose and boundary

The model exchange request intentionally contains stable IDs and SHA-256
fingerprints instead of embedding HostState or PhageProfile records. This
registry resolves those exact references to local JSON records.

The registry is a local lookup and integrity mechanism. An entry marked
`available` means only that the local operator has registered the file for
resolution. It does **not** certify biological validity, source quality, assay
compatibility, sequence correctness, or fitness for model training.

This milestone resolves HostState and PhageProfile records only. It does not
download remote content, resolve raw host/phage genome files, derive features,
train models, or run inference. The returned summary explicitly records that
sequence assets remain unresolved and biological validity has not been
established.

## Manifest format

A registry is a versioned JSON document with exactly these top-level fields:

- `schema_version`: currently `1`;
- `registry_id`: a non-empty local registry identifier;
- `entries`: an array of exact resource mappings.

The machine-readable schema is
[`resource-registry-v1.schema.json`](../schemas/resource-registry-v1.schema.json).

Each entry contains exactly:

- `kind`: `host_state` or `phage_profile`;
- `resource_id`: HostState `record_id` or PhageProfile `profile_id`;
- `content_sha256`: canonical record fingerprint used by ModelRequest;
- `file_sha256`: SHA-256 of the exact bytes of the JSON file on disk;
- `path`: a POSIX-style path relative to the registry file's directory;
- `status`: `available` or `withdrawn`.

An identifier may appear only once per resource kind, and a file path may appear
only once. Resource paths cannot be absolute, contain `.` or `..` components,
use backslashes, or resolve through a symlink outside the registry directory.
The registry never follows a path supplied by ModelRequest; all paths come from
the locally managed manifest.

## Integrity checks

For each resource, the resolver:

1. validates ModelRequest v1 and its canonical content hash;
2. finds the exact kind and resource ID in the registry;
3. requires `status = available` and an exact match for the request's content
   SHA-256;
4. confines the resolved path to the registry directory;
5. checks the file's raw-byte SHA-256 against `file_sha256`;
6. parses a JSON object and recomputes its canonical content SHA-256, excluding
   an optional embedded `content_sha256` field;
7. checks the resource's identity fields against ModelRequest.

HostState checks include record ID, HostProfile ID, genome fingerprint and
`as_of`. PhageProfile checks include profile ID, revision, schema version and
its canonical fingerprint. Missing, withdrawn, tampered, ambiguous, or mismatched
resources fail closed.

The manifest itself is also hashed by exact file bytes and that digest is
included in the resolver summary. This is useful for audit and reproducibility;
it is not a digital signature and does not prove who created or approved the
manifest. Keep production registries under controlled versioning and review.

## Usage

Create a registry only after the exact HostState and PhageProfile JSON records
have been exported and reviewed. Put their files under the registry directory,
compute both hashes, and add one entry for each record. No production registry
or biological record is bundled with this repository.

Then run:

```bash
phage-model-interface validate-request path/to/model-request.json

phage-model-interface resolve-request-records \
  path/to/model-request.json \
  --registry path/to/registry.v1.json
```

The command prints a compact audit summary containing the request hash, registry
hash, resolved record IDs, content/file hashes and relative paths. It does not
print or copy the resource records themselves.

## Next boundary

The sequence-resource layer is specified in [sequence-registry.md](sequence-registry.md), and the downstream provenance manifest is defined in [FeatureArtifact v1](feature-artifact.md). Successful record or sequence resolution does not run an analyzer, derive features, or establish biological validity.
