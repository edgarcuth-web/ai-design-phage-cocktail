"""Download and verify the pinned GenoPHI E. coli subset."""

from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen

from phage_cocktail.registry import load_dataset_manifest

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPOSITORY_ROOT / "data/registry/genophi_ecoli_subset.v1.json"


def main() -> None:
    manifest = load_dataset_manifest(MANIFEST_PATH)
    source_url = manifest["source"]["download_url"]
    output_path = REPOSITORY_ROOT / manifest["local_path"]

    with urlopen(source_url, timeout=60) as response:
        payload = response.read()

    observed_sha256 = sha256(payload).hexdigest()
    if observed_sha256 != manifest["sha256"]:
        raise SystemExit(
            "Downloaded dataset checksum mismatch. "
            f"Expected {manifest['sha256']}, observed {observed_sha256}."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_path.write_bytes(payload)
    temporary_path.replace(output_path)

    print(f"Downloaded {len(payload):,} bytes")
    print(f"SHA-256: {observed_sha256}")
    print(f"Source: {source_url}")
    print(f"Saved to: {output_path}")


if __name__ == "__main__":
    main()
