import pytest

from phage_cocktail.baseline import cocktail_success_probability


def test_single_phage_probability_is_unchanged():
    assert cocktail_success_probability([0.7]) == pytest.approx(0.7)


def test_independent_cocktail_probability():
    assert cocktail_success_probability([0.5, 0.5]) == pytest.approx(0.75)


def test_certain_success_if_any_probability_is_one():
    assert cocktail_success_probability([0.2, 1.0, 0.1]) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "values",
    [
        [],
        [-0.1],
        [1.1],
        [float("nan")],
        [float("inf")],
    ],
)
def test_invalid_probabilities_raise(values):
    with pytest.raises(ValueError):
        cocktail_success_probability(values)
