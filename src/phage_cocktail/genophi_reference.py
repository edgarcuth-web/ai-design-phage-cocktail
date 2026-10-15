"""Helpers for running GenoPHI as an external reference baseline."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

GENOPHI_VERSION = "1.0.3"
GENOPHI_REPOSITORY = "https://github.com/Noonanav/GenoPHI.git"
GENOPHI_COMMIT = "f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f"


def build_demo_command(
    data_root: str | Path,
    output_dir: str | Path,
    *,
    threads: int = 4,
) -> list[str]:
    """Build the authors' small protein-family demo command.

    The demo contains 25 E. coli strains, 25 phages, and 625 measured
    interactions. filter_type=strain keeps strains separated during
    GenoPHI's repeated train/test evaluation.
    """
    data_root = Path(data_root)
    output_dir = Path(output_dir)

    return [
        "genophi",
        "protein-family-workflow",
        "--input_path_strain",
        str(data_root / "strain_AAs"),
        "--input_path_phage",
        str(data_root / "phage_AAs"),
        "--phenotype_matrix",
        str(data_root / "ecoli_test_interaction_matrix.csv"),
        "--output_dir",
        str(output_dir),
        "--threads",
        str(threads),
        "--num_features",
        "50",
        "--num_runs_fs",
        "5",
        "--num_runs_modeling",
        "10",
        "--method",
        "rfe",
        "--filter_type",
        "strain",
    ]


def find_performance_file(output_dir: str | Path) -> Path:
    """Locate GenoPHI's summary model-performance CSV."""
    output_dir = Path(output_dir)
    candidates = [
        output_dir
        / "modeling_results"
        / "model_performance"
        / "model_performance_metrics.csv",
        output_dir
        / "modeling_results"
        / "select_features_model_performance.csv",
    ]

    for path in candidates:
        if path.exists():
            return path

    tried = "\n".join(str(path) for path in candidates)
    raise FileNotFoundError(f"no GenoPHI performance file found; tried:\n{tried}")


def summarize_performance(path: str | Path) -> pd.DataFrame:
    """Return mean/std for useful numeric metrics from a GenoPHI result CSV."""
    frame = pd.read_csv(path)
    numeric = frame.select_dtypes(include="number")

    preferred = [
        column
        for column in numeric.columns
        if any(
            token in column.lower()
            for token in (
                "auc",
                "precision",
                "recall",
                "f1",
                "mcc",
                "accuracy",
                "average_precision",
            )
        )
    ]

    selected = numeric[preferred] if preferred else numeric
    if selected.empty:
        raise ValueError("GenoPHI performance CSV contains no numeric metrics")

    return selected.agg(["mean", "std"]).T
