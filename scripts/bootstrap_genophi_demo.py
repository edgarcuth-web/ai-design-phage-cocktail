"""Fetch only GenoPHI's pinned test-data directory using sparse checkout."""

from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

from phage_cocktail.genophi_reference import GENOPHI_COMMIT, GENOPHI_REPOSITORY

TARGET = Path("external/GenoPHI")


def run(*args: str) -> None:
    subprocess.run(args, check=True)


def current_commit() -> str | None:
    if not (TARGET / ".git").exists():
        return None

    result = subprocess.run(
        ["git", "-C", str(TARGET), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def main() -> None:
    if shutil.which("git") is None:
        raise SystemExit("git is required but was not found on PATH")

    if TARGET.exists() and not (TARGET / ".git").exists():
        raise SystemExit(
            f"{TARGET} already exists but is not a Git checkout; move or remove it first"
        )

    if not TARGET.exists():
        TARGET.parent.mkdir(parents=True, exist_ok=True)
        run(
            "git",
            "clone",
            "--filter=blob:none",
            "--no-checkout",
            GENOPHI_REPOSITORY,
            str(TARGET),
        )
        run("git", "-C", str(TARGET), "sparse-checkout", "init", "--cone")
        run(
            "git",
            "-C",
            str(TARGET),
            "sparse-checkout",
            "set",
            "data/test_data",
        )

    if current_commit() != GENOPHI_COMMIT:
        run(
            "git",
            "-C",
            str(TARGET),
            "fetch",
            "--depth",
            "1",
            "origin",
            GENOPHI_COMMIT,
        )
        run(
            "git",
            "-C",
            str(TARGET),
            "checkout",
            "--detach",
            "FETCH_HEAD",
        )

    data_root = TARGET / "data" / "test_data"
    required = [
        data_root / "strain_AAs",
        data_root / "phage_AAs",
        data_root / "ecoli_test_interaction_matrix.csv",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        names = "\n".join(str(path) for path in missing)
        raise SystemExit(f"GenoPHI checkout is missing expected demo data:\n{names}")

    print(f"GenoPHI demo data ready at: {data_root}")
    print(f"Pinned commit: {GENOPHI_COMMIT}")


if __name__ == "__main__":
    main()
