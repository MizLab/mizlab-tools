from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from mizlab_tools import gbk_utils


def load(gbkname) -> Path:
    """test_get_creature_name.

    Args:
        gbk:
        expected:
    """
    return Path(Path(__file__).resolve().parent / "in" / gbkname)


def make_mock_record(seq: Seq = Seq("ATGC"),
                     id: str = "mock_id",
                     name: Optional[str] = "mock_name",
                     description: Optional[str] = "mock_description",
                     dbxrefs: Optional[List[str]] = ["test_dbxref"],
                     features: Optional[List[SeqFeature]] = [SeqFeature()],
                     annotations: Dict[str, str] = {"": ""}) -> SeqRecord:
    """make_mock_record.

    Args:
        seq (Seq): seq
        id (str): id
        name (Optional[str]): name
        description (Optional[str]): description
        dbxrefs (Optional[List[str]]): dbxrefs
        features (Optional[List[SeqFeature]]): features
        annotations (Dict[str, str]): annotations

    Returns:
        SeqRecord:
    """
    return SeqRecord(seq=seq,
                     id=id,
                     name=name,
                     description=description,
                     dbxrefs=dbxrefs,
                     features=features,
                     annotations=annotations)


@pytest.mark.parametrize(("qualifiers_key", "db_xref", "expected"),
                         [("db_xref", "taxon:45628", "45628"),
                          ("db_xref", "Taxon:45628:hogehoge", "45628"),
                          ("not_db_xref", "taxon:1111", None),
                          ("db_xref", "not_taxon:1111", None),
                          ("not_db_xref", "not_taxon:1111", None)])
def test_get_taxonID(qualifiers_key, db_xref, expected):
    mock_feature = SeqFeature(FeatureLocation(1, 100),
                              type="source",
                              qualifiers={f"{qualifiers_key}": [f"{db_xref}"]})
    mock_record = make_mock_record(features=[mock_feature])
    if expected is not None:
        assert gbk_utils.get_taxonID(mock_record) == expected
    else:
        assert gbk_utils.get_taxonID(mock_record) is None


@pytest.mark.parametrize(("source", "expected"),
                         [("mock_definition", "mock_definition")])
def test_get_definition(source: str, expected: str):
    assert gbk_utils.get_definition(make_mock_record(description=source))


@pytest.mark.parametrize(("k", "v", "expected"), [("organism", "test", "test"),
                                                  ("organisms", "test", None)])
def test_get_creature_name(k: str, v: str, expected: Optional[str]):
    if expected is None:
        assert gbk_utils.get_creature_name(make_mock_record(annotations={k: v})) is None
    else:
        assert gbk_utils.get_creature_name(make_mock_record(annotations={k: v})) == v


@pytest.mark.parametrize(("source", "window_size", "overhang"), [
    ("abcde", 3, None),
    (Seq("ATGCATGC"), 4, None),
    ("abcde", 3, "before"),
    ("abcde", 3, "after"),
    ("abcde", 3, "both"),
    ("abcde", 6, None),
])
def test_window_search(source: Union[str, Seq], window_size: int,
                       overhang: Optional[gbk_utils.OVERHANG]):
    if overhang is None:
        if len(source) < window_size:
            for w in gbk_utils.window_search(source, window_size, overhang):
                assert False    # this line must not go througth
            else:
                assert True
        else:
            assert all([
                lambda x: len(x) == window_size
                for w in gbk_utils.window_search(source, window_size, overhang)
            ])
    else:
        windows = [w for w in gbk_utils.window_search(source, window_size, overhang)]
        if overhang in {"before", "both"}:
            for i in range(window_size - 1):
                assert len(windows[i]) == i + 1
        if overhang in {"after", "both"}:
            for i in range(1, window_size):
                assert len(windows[-i]) == i


@pytest.mark.parametrize(("name", "expected"), [
    (None, False),
    ("", False),
    ("foo x bar", True),
    ("foo  x bar", True),
    ("foo x  bar", True),
    ("foox bar", False),
    ("foo xbar", False),
    ("fooxbar", False),
])
def test_is_mongrel(name: Optional[str], expected: bool):
    if name is None:
        mock = make_mock_record()
    else:
        mock = make_mock_record(annotations={"organism": name})
    assert gbk_utils.is_mongrel(mock) == expected


@pytest.mark.parametrize(
    ("source", "expected"),
    [("join(LVXP01042324.1:1..16300)", ("LVXP01042324", 0, 16300, False)),
     ("join(JABSTT010003572.1:1..14723)", ("JABSTT010003572", 0, 14723, False)),
     ("join(complement(JAACYO010019948.1:1..16778))",
      ("JAACYO010019948", 0, 16778, True)), ("this is not contig", None)],
)
def test_parse_contig(source: str, expected: Optional[Union[Tuple[str, int, int,
                                                                  bool]]]):
    if expected is None:
        assert gbk_utils.parse_contig(source) is None
    else:
        key = ("accession", "start", "end", "is_complement")
        expected_dict = {k: v for k, v in zip(key, expected)}
        assert gbk_utils.parse_contig(source) == expected_dict


#
# @pytest.mark.parametrize(("gbk", "expected"), [("NC_012920.gbk", True),
#                                                ("no_seq_creature.gbk", False)])
# def test_has_seq(gbk, expected):
#     assert gbk_utils.has_seq(load(gbk)) == expected
#
#
# @pytest.mark.parametrize(("name", "expected"), [("Homo sapiens", False),
#                                                 ("foo x bar", True),
#                                                 ("includexinname", False),
#                                                 ("include xinname", False)])
# def test_is_mongrel(name, expected):
#     assert gbk_utils.is_mongrel(name) == expected
#
#
# @pytest.mark.parametrize(("definition", "expected"), [
#     ("Homo sapiens mitochondrion, complete genome.", True),
#     ("Homo sapiens isolate CHM1 mitochondrion, complete sequence, whole genome shotgun sequence.",
#      False),
#     ("Kudoa iwatai mitochondrion, chromosome 2, complete sequence.", False),
# ])
# def test_is_complete_genome(definition, expected):
#     assert gbk_utils.is_complete_genome(definition) == expected
#
#
# @pytest.mark.parametrize(("source", "expected"),
#                          [("ATGCATGC", "ATGCATGC"), ("AAAANNNN", "AAAA"),
#                           ("aTaT", "aTaT"), ("aaaannnn", "aaaa")])
# def test_to_only_atgc(source, expected):
#     assert gbk_utils.to_only_actg(source) == Seq(expected)
#
#
# @pytest.mark.parametrize(
#     ("source", "expected"),
#     [("join(LVXP01042324.1:1..16300)", ("LVXP01042324", 0, 16300, False)),
#      ("join(JABSTT010003572.1:1..14723)", ("JABSTT010003572", 0, 14723, False)),
#      ("join(complement(JAACYO010019948.1:1..16778))",
#       ("JAACYO010019948", 0, 16778, True))])
# def test_parse_contig(source, expected):
#     assert gbk_utils.parse_contig(source) == {
#         k: v
#         for k, v in zip(("accession", "start", "end", "is_complement"), expected)
#     }
#
#
# @pytest.mark.parametrize(("source", "expecteds"), [("NC_046603.gbk", (True, )),
#                                                    ("NC_012920.gbk", (False, ))])
# def test_has_contig(source, expecteds):
#     for record, expected in zip(SeqIO.parse(load(source), "genbank"), expecteds):
#         assert gbk_utils.has_contig(record) == expected
#
#
# @pytest.mark.parametrize(("source", "expecteds"), [("NC_046603.gbk", (False, )),
#                                                    ("NC_012920.gbk", (True, ))])
# def tese_get_seq(source, expecteds):
#     for record, expected in zip(SeqIO.parse(load(source), "genbank"), expecteds):
#         assert (record.seq == gbk_utils.get_seq(record,
#                                                 recursive=True,
#                                                 search_gbk_root=load(""))) == expected
