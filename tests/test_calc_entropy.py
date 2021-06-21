import pytest

from mizlab_tools import calculate_weights


@pytest.mark.parametrize(("source", "expected"), [(("ATGC", ), ({"ATG": 1, "TGC": 1}))])
def test_count_triplets(source, expected):
    assert calculate_weights.count_words(source, "ATGC", 3) == expected


@pytest.mark.parametrize(("source", "expected"), [({
    "AAA": 1,
    "AAT": 1,
    "AAG": 1,
    "AAC": 1
}, {
    "AAA": 2,
    "AAT": 2,
    "AAG": 2,
    "AAC": 2
})])
def test_calc_self_entropy(source, expected):
    assert calculate_weights.compute_self_entropy(source, "ATGC") == expected
