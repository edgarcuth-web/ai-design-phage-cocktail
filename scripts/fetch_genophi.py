"""Download the pinned GenoPHI E. coli subset used for the first baseline."""

from pathlib import Path
from urllib.request import urlopen

GENOPHI_COMMIT = "f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f"
SOURCE_URL = (
    "https://raw.githubusercontent.com/Noonanav/GenoPHI/"
    f"{GENOPHI_COMMIT}/data/interaction_matrices/"
    "ecoli_interaction_matrix_subset.csv"
)
OUTPUT_PATH = Path("data/raw/genophi_ecoli_subset.csv")


def main() -> None:
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    with urlopen(SOURCE_URL, timeout=60) as response:
        payload = response.read()

    OUTPUT_PATH.write_bytes(payload)
    print(f"Downloaded {len(payload):,} bytes")
    print(f"Source: {SOURCE_URL}")
    print(f"Saved to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
