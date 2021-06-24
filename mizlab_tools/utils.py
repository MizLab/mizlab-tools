from itertools import groupby
from typing import Generator, Iterable, Iterator, TypeVar

T = TypeVar("T")


def split_per_n(it: Iterable[T], n: int) -> Iterator[Generator[T, None, None]]:
    """Split iterable object per N.

    Args:
        it (Iterable[T]): Iterable object.
        n (int): length of splited size.

    Returns:
        Iterator[Generator[T, None, None]]: Splited object.
    """
    for _, item in groupby(enumerate(it), lambda x: x[0] // n):
        yield (x[1] for x in item)
