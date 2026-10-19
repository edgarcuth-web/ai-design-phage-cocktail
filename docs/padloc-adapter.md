# PADLOC host-defense adapter v1

## Purpose and scope

The first real feature analyzer is [PADLOC](https://github.com/padlocbio/padloc), a tool that searches prokaryotic genomes against profile HMMs and system-architecture models to identify known antiviral defense systems. This adapter invokes the external PADLOC executable on a HostState's exact, registry-verified host-genome FASTA.

This is a real-tool adapter, not a simulation. It does not contain fake biological calls and does not manufacture a positive result when PADLOC reports no systems. It does not automatically promote output rows into HostProfile components or make phage-susceptibility predictions.

The adapter intentionally pins both PADLOC and PADLOC-DB to `v2.0.0`. PADLOC's maintainers explicitly state that PADLOC 2.x and PADLOC-DB 2.x must be used together; the exact versions are checked before execution. The installed database tree is also content-addressed, so two directories bearing the same release label but different files do not share a digest.

## Installation

Use the committed Conda environment:

```bash
conda env create -f envs/padloc.yml
conda activate phage-padloc
python -m pip install -e .
```

Install the exact database release into a dedicated path that contains no unrelated files:

```bash
mkdir -p "$HOME/.local/share/phage-cocktail"
PADLOC_DATA="$HOME/.local/share/phage-cocktail/padloc-db-v2.0.0"
padloc --data "$PADLOC_DATA" --db-install v2.0.0
padloc --version
padloc --data "$PADLOC_DATA" --db-version
padloc --check-deps
```

**Important:** upstream PADLOC's `--db-install` removes the directory passed through `--data` if that directory already exists, then replaces it with the database. Use only a dedicated database directory. The adapter itself never downloads, updates, compiles, or modifies PADLOC-DB.

The expected version output is `padloc v2.0.0` and `padloc-db v2.0.0`. Missing runtime dependencies are treated as an error, not silently ignored.

## Input contract

The command requires:

- A valid HostState v1 JSON record with a correct canonical `content_sha256` and a host-genome SHA-256.
- A local sequence registry containing exactly one available `host_genome` entry matching that normalized sequence digest. The adapter verifies the registry, exact FASTA file-byte digest, normalized sequence digest, FASTA structure, base count and contig count.
- A PADLOC data directory with the exact v2.0.0 release marker and required compiled database/model files.

A PhageProfile and pairwise ModelRequest are deliberately **not** required: bacterial defense annotation is a host-level feature step and should be reusable independently of any candidate phage. The HostState ID, HostState content digest, HostState raw-file digest, sequence-registry ID/digest, and both sequence fingerprints are carried into the result provenance.

The FASTA normalization convention is defined by [the sequence registry contract](sequence-registry.md). A normalized sequence digest identifies the sequence content independently of FASTA header text and wrapping; the raw-file digest still records the exact bytes consumed from the registry.

## Run

```bash
phage-model-interface annotate-host-defense-padloc \
  --host-state path/to/host-state.v1.json \
  --sequence-registry path/to/sequence-registry.v1.json \
  --padloc-data "$PADLOC_DATA" \
  --outdir runs/padloc-host-001 \
  --cpu 4
```

The output directory must be new or empty. This prevents an old `*_padloc.csv` from being mistaken for output from the current run. The adapter stages the verified FASTA under a temporary, whitespace-free path because the upstream PADLOC shell wrapper does not quote one internal basename expansion. It does not alter the registered input FASTA.

The adapter executes the native `padloc --fna ... --outdir ... --cpu ...` path. PADLOC runs its own Prodigal gene-calling workflow before HMMER/system-model analysis. The adapter does not supply separate CRISPR-array or ncRNA input files in this first integration.

## Outputs

A successful run creates:

- `feature-artifact.v1.json`: the general FeatureArtifact manifest, including input/record lineage, exact PADLOC executable digest, canonical parameter digest, timestamps, pinned database version and tree digest, output digests, dependency-check evidence and limitations.
- `outputs/host-defense-summary.json`: the versioned PADLOC Host Defense Result v1 summary, containing the native CSV rows when present, distinct system count, gene-hit count, system names and the native output path.
- `outputs/padloc-native/`: native PADLOC outputs, including the CSV/GFF when systems are found and intermediate Prodigal/HMMER files produced by the tool.
- `logs/`: PADLOC version, database version, dependency check, standard output and standard error.

The database tree digest is calculated as SHA-256 over canonical JSON containing the sorted relative paths, SHA-256 digests and byte sizes of all regular files in the installed database directory. Its method identifier is stored in the result. This is a digest of the **installed/compiled tree**, not a claim that the upstream release archive itself was hashed.

Verify the manifest and every referenced file:

```bash
phage-model-interface verify-feature-artifact runs/padloc-host-001/feature-artifact.v1.json
```

Validate the structured PADLOC summary in Python:

```python
from phage_cocktail.padloc_adapter import load_padloc_summary

summary = load_padloc_summary(
    "runs/padloc-host-001/outputs/host-defense-summary.json"
)
print(summary["outcome"], summary["systems_detected"], summary["gene_hit_count"])
```

## Interpretation and limitations

`hits_detected` means PADLOC emitted one or more rows in its native result. `no_hits_detected` means PADLOC completed successfully and emitted its explicit no-hit completion message without a result CSV. A missing CSV without that completion message is treated as an error. A zero-row CSV without the explicit no-hit message is also treated as an error.

**No hits do not prove that the bacterium has no anti-phage defenses.** This workflow is limited to systems represented in PADLOC-DB v2.0.0, the genome/assembly supplied, PADLOC's detection thresholds, and the `--fna` path's gene predictions. Some systems may be novel, divergent, incomplete, or depend on RNA/CRISPR evidence not supplied here. A detected system is a computational annotation, not an experimentally demonstrated phenotype.

Do not map a missing row to a biological negative. Before using these calls in HostProfile or CocktailDesign, the next layer should preserve the exact source row, gene/contig coordinates, HMM identity, E-values/coverage fields, database release and the distinction between “not detected,” “not searched,” and “unknown.” Orthogonal evidence and experimental phenotypes remain separate.

## Reproducibility and tests

The unit tests use a controlled executable stub to exercise version checks, error handling, CSV parsing, zero-hit semantics and provenance. The PADLOC integration workflow additionally installs PADLOC v2.0.0 and PADLOC-DB v2.0.0, downloads a test genome from a pinned commit of the upstream PADLOC repository, runs the real analyzer, and verifies the resulting FeatureArtifact. The upstream test genome and outputs are not committed to this repository.

## Citation and license

Please cite:

Payne, L. J., Todeschini, T. C., Wu, Y., Perry, B. J., Ronson, C. W., Fineran, P. C., Nobrega, F. L., Jackson, S. A. (2021). Identification and classification of antiviral defence systems in bacteria and archaea with PADLOC reveals new system types. *Nucleic Acids Research*, 49, 10868–10878. https://doi.org/10.1093/nar/gkab883

PADLOC and PADLOC-DB repositories identify their software/database license as MIT; retain upstream copyright and license notices when redistributing their files.
