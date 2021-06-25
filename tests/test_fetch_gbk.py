from typing import Iterable

import pytest

from mizlab_tools import fetch_gbk


@pytest.mark.parametrize(("accessions"), [
    (["NC_012920", "NC_005089"]),
])
def test_fetch(accessions: Iterable[str]):
    for accession, r in zip(accessions,
                            fetch_gbk.fetch(accessions, email="example@example.com")):
        assert r.name == accession
