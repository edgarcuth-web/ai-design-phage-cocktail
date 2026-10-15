import pandas as pd

from phage_cocktail.analysis import summarize_interactions


def test_summarize_complete_interaction_matrix():
    frame = pd.DataFrame(
        {
            "bacterium": ["B1", "B1", "B2", "B2"],
            "phage": ["P1", "P2", "P1", "P2"],
            "infection": [1, 0, 0, 1],
        }
    )

    summary = summarize_interactions(frame)

    assert summary["interactions"] == 4
    assert summary["bacteria"] == 2
    assert summary["phages"] == 2
    assert summary["positive"] == 2
    assert summary["negative"] == 2
    assert summary["positive_rate"] == 0.5
    assert summary["observed_fraction"] == 1.0
