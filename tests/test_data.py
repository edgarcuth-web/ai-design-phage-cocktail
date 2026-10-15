import pandas as pd
import pytest

from phage_cocktail.data import load_interactions


def test_load_interactions_accepts_binary_table(tmp_path):
    path = tmp_path / "interactions.csv"
    pd.DataFrame(
        {
            "bacterium": ["B1", "B1", "B2"],
            "phage": ["P1", "P2", "P1"],
            "infection": [1, 0, 0],
        }
    ).to_csv(path, index=False)

    frame = load_interactions(path)

    assert list(frame.columns) == ["bacterium", "phage", "infection"]
    assert len(frame) == 3


def test_load_interactions_rejects_nonbinary_labels(tmp_path):
    path = tmp_path / "interactions.csv"
    pd.DataFrame(
        {
            "bacterium": ["B1"],
            "phage": ["P1"],
            "infection": [2],
        }
    ).to_csv(path, index=False)

    with pytest.raises(ValueError, match="binary"):
        load_interactions(path)
