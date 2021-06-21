import os
import re
from collections import Counter
from pathlib import Path
from typing import (Any, AnyStr, Dict, Iterable, Iterator, Optional, Sequence,
                    Tuple, TypeVar, Union)

from Bio import Seq, SeqIO, SeqRecord

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

T = TypeVar("T")
Openable = Union[str, bytes, int, "os.PathLike[Any]"]


def get_taxonID(path: Openable) -> str:
    """対象のgbkファイルのrecordからtaxonIDを抽出する

    Args:
        record (Openable): 対象のファイルのpath

    Returns:
        str: db_xrefに記載されたTaxonID
    """
    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                db_xref: str = feature.qualifiers["db_xref"][0]
                taxonID: str = db_xref.split(":")[1]
                return taxonID
        raise NotFoundTaxonIDError(f"Can not Found taxonID in {record.name}")
    return ""


class NotFoundTaxonIDError(Exception):
    pass


def get_definition(path: Openable) -> str:
    """get definition.

    Args:
        path (Path): 対象のファイルのpath

    Returns:
        str: 生物の学名
    """
    for record in SeqIO.parse(path, "genbank"):
        return record.description
    return ""


def get_creature_name(path: Openable) -> str:
    """get creature name.

    Args:
        path (Openable): 対象のファイルのpath

    Returns:
        str:
    """
    for record in SeqIO.parse(path, "genbank"):
        return record.annotations["organism"]
    return ""


OVERHANG = Literal["before", "after", "both"]


def window_search(target: Union[str, Seq.Seq],
                  window_size: int,
                  overhang: Optional[OVERHANG] = None) -> Iterator[Union[str, Seq.Seq]]:
    """window_search.

    Args:
        target (Iterable): list like object.
        window_size (int): window_size
        overhang (Optional[OVERHANG]): overhang

    Returns:
        Iterator[T]:
    """

    if overhang in {"before", "both"}:
        for i in range(1, window_size):
            yield target[:i]
    for i in range(len(target) - window_size + 1):
        yield target[i:i + window_size]
    if overhang in {"after", "both"}:
        for i in range(len(target) - window_size + 1, len(target)):
            yield target[i:]


# def get_rate(string: AnyStr, allowed: AnyStr) -> Dict[str, float]:
#     pretty = re.sub(f"[^{allowed}]", "", string.upper())
#     length = len(pretty)
#     counter = Counter(pretty)
#     rate = {k: v / lengtg for k, v in counter.items()}
#     return rate


def has_seq(gbk: Openable) -> bool:
    """与えられたgbkファイルが有効な配列長を持つかどうかを返す.

    Args:
        gbk (Openable): gbkファイルへのpath

    Returns:
        bool:
    """
    return any([
        len((re.sub("[^ATGC]", "",
                    str(rec.seq).upper()))) for rec in SeqIO.parse(gbk, "genbank")
    ])


def is_mongrel(name: str) -> bool:
    """'~ x ~'で書かれる雑種かどうかを返す.

    Args:
        name (str): 生物種

    Returns:
        bool:
    """
    return " x " in name


def is_complete_genome(definition: str) -> bool:
    """完全なミトコンドリアゲノムかどうかを返す.

    Args:
        definition (str): definition

    Returns:
        bool:
    """
    return "complete" in definition and "shotgun" not in definition and "chromosome" not in definition


def parse_contig(contig: str) -> Optional[Dict[str, Union[str, int]]]:
    """return contig doscription

    Args:
        contig [str]: A contig string.

    Returns:
        dict: A dict of contig informations.
    """
    contig_pattern = re.compile(
        r"join(\(complement)?\((\w*(\.\d)?):(\d+)\.\.(\d+)\)\)?")
    match_result = contig_pattern.match(contig)
    if match_result is None:
        return None
    else:
        group = match_result.groups()
        is_complement = group[0] is not None
        accession = group[1].split(".")[0]
        start = int(group[3]) - 1
        end = int(group[4])
        return {
            "accession": accession,
            "is_complement": is_complement,
            "start": start,
            "end": end
        }


def has_contig(record: SeqRecord) -> bool:
    """Does the record has contig?

    Args:
        record (SeqRecord): record 

    Returns:
        bool: have contig or not have.
    """
    return "contig" in record.annotations


def get_seq(record: SeqRecord,
            recursive: bool = False,
            search_gbk_root: Optional[Openable] = None,
            is_complement: bool = False) -> Seq.Seq:
    """gbkから配列を取得する
    recursiveがTrueの時、contigに書かれたものを取得する

    Args:
        record (SeqRecord): 取得対象のレコード
        recursive (bool, optional): contigがあったときに再帰的に取得するか. Defaults to False.
        search_gbk_root (Optional[Openable], optional): recursiveがTrueの時、どこにあるgbkを探す対象にするか. Defaults to None.
        is_complement (bool): 再帰的に見る時

    Raises:
        FileNotFoundError: recursiveがtrueで探したが、contigのgbkが見つからなかったときに送出される。

    Returns:
        Seq.Seq: record's sequence.
    """
    if recursive and has_contig(record):
        contig_info = parse_contig(record.annotations["contig"])
        contig_gbk = Path(search_gbk_root) / f"{contig_info['accession']}.gbk"
        if not contig_gbk.exists():
            raise FileNotFoundError
        else:
            for r in SeqIO.parse(contig_gbk, "genbank"):    # 複数レコードは考慮しない(面倒なので)
                return get_seq(r,
                               recursive=True,
                               search_gbk_root=search_gbk_root,
                               is_complement=(is_complement
                                              ^ contig_info["is_complement"]))
                # よって複数レコードだと最初のものが対象になるが多分大丈夫
    else:
        seq = record.seq
        if is_complement:
            return seq.complement()
        else:
            return seq


class FileNotFoundError(Exception):
    pass
