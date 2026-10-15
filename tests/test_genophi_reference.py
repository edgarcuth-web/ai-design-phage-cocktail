import pandas as pd
import pytest

from phage_cocktail.genophi_reference import (
    build_demo_command,
    find_performance_file,
    summarize_performance,
)


def test_build_demo_command_enforces_strain_filtered_evaluation():
    command = build_demo_command("demo_data", "results", threads=6)

    assert command[0:2] == ["genophi", "protein-family-workflow"]
    assert "--filter_type" in command
    assert command[command.index("--filter_type") + 1] == "strain"
    assert command[command.index("--threads") + 1] == "6"
    assert command[command.index("--num_features") + 1] == "50"


def test_find_and_summarize_performance(tmp_path):
    performance_dir = tmp_path / "modeling_results" / "model_performance"
    performance_dir.mkdir(parents=True)
    path = performance_dir / "model_performance_metrics.csv"

    pd.DataFrame(
        {
            "AUC": [0.7, 0.9],
            "Accuracy": [0.6, 0.8],
            "run": [1, 2],
        }
    ).to_csv(path, index=False)

    found = find_performance_file(tmp_path)
    summary = summarize_performance(found)

    assert found == path
    assert summary.loc["AUC", "mean"] == pytest.approx(0.8)
    assert summary.loc["Accuracy", "mean"] == pytest.approx(0.7)
    assert "run" not in summary.index
