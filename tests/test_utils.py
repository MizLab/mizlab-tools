from typing import Any, Iterable
from mizlab_tools import utils
import pytest


@pytest.mark.parametrize(("source", "org_len", "n"), [
    ("hogehgoe", 8, 3),
    ([1, 2, 3, 4, 5], 5, 2),
    ((i for i in range(10)), 10, 4),
])
def test_split_per_n(source: Iterable[Any], org_len: int, n: int):
    for splitted in utils.split_per_n(source, n):
        fixed = tuple(splitted)
        assert len(fixed) == n or len(fixed) == org_len % n
