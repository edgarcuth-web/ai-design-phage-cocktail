"""Run GenoPHI's official 25 x 25 protein-family demo."""

from __future__ import annotations

from pathlib import Path
import shlex
import shutil
import subprocess

from phage_cocktail.genophi_reference import build_demo_command

DATA_ROOT = Path("external/GenoPHI/data/test_data")
OUTPUT_DIR = Path("artifacts/genophi_demo")


def main() -> None:
    if not DATA_ROOT.exists():
        raise SystemExit(
            "GenoPHI demo data is missing. Run: "
            "python scripts/bootstrap_genophi_demo.py"
        )

    if shutil.which("genophi") is None:
        raise SystemExit(
            "The genophi executable was not found. "
            'Install the reference extra with: pip install -e ".[reference]"'
        )

    if shutil.which("mmseqs") is None:
        raise SystemExit(
            "MMseqs2 was not found on PATH. Install MMseqs2 before running "
            "the protein-family reference workflow."
        )

    OUTPUT_DIR.parent.mkdir(parents=True, exist_ok=True)
    command = build_demo_command(DATA_ROOT, OUTPUT_DIR, threads=4)

    print("Running:")
    print(shlex.join(command))
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
