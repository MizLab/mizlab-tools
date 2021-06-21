from itertools import groupby
from typing import Iterable, Iterator, TypeVar

T = TypeVar("T")


def split_per_n(it: Iterable[T], n: int) -> Iterator[T]:
    for _, item in groupby(enumerate(it), lambda x: x[0] // n):
        yield (x[1] for x in item)
