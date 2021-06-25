from typing import Dict, List, Optional, Union

import pytest

from mizlab_tools import calculate_coordinates
from tests.test_gbk_utils import make_mock_record


@pytest.mark.parametrize(("seq", "mapping", "weight", "expected"), [
    ("ATGC", {
        "A": [1, 1],
        "T": [-1, 1],
        "G": [-1, -1],
        "C": [1, -1]
    }, None, [[0, 0], [1, 1], [0, 2], [-1, 1], [0, 0]]),
    ("ATGC", {
        "A": [0],
        "T": [0],
        "G": [0],
        "C": [0]
    }, None, [[0] for i in range(5)]),
    ("ATGC", {
        "A": [1],
        "T": [1],
        "G": [1],
        "C": [1]
    }, {
        "ATG": 2,
        "TGC": 0
    }, [[0], [1], [2], [4], [4]]),
])
def test_calc_coord(
    seq: str,
    mapping: Dict[str, Union[int, float]],
    weight: Optional[Dict[str, Union[int, float]]],
    expected: List[List[Union[int, float]]],
):
    assert calculate_coordinates.calc_coord(make_mock_record(seq=seq), mapping,
                                            weight) == expected
