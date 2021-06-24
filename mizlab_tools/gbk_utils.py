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


def get_taxonID(record: SeqRecord.SeqRecord) -> Optional[str]:
    """get taxonID.

    Args:
        record (SeqRecord.SeqRecord): seq record.

    Returns:
        Optional[str]: if exists taxon id. else None.
    """
    for feature in record.features:
        if feature.type == "source":
            if feature.qualifiers.get("db_xref", None) is not None:
                for db_xref in feature.qualifiers["db_xref"]:
                    name, identifier, *_ = db_xref.split(":")
                    if name.lower() == "taxon":
                        return identifier
    return None


def get_definition(record: SeqRecord.SeqRecord) -> str:
    """get definition.

    Args:
        record (SeqRecord.SeqRecord): record

    Returns:
        str: definition string.
    """
    return record.description


def get_creature_name(record: SeqRecord.SeqRecord) -> Optional[str]:
    """get creature name.

    Args:
        record (SeqRecord.SeqRecord): record

    Returns:
        Optional[str]:
    """
    return record.annotations.get("organism", None)


OVERHANG = Literal["before", "after", "both"]


def window_search(target: Union[str, Seq.Seq],
                  window_size: int,
                  overhang: Optional[OVERHANG] = None) -> Iterator[Union[str, Seq.Seq]]:
    """window_search.

    Args:
        target (Iterable): like sequence.
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


def has_seq(record: SeqRecord.SeqRecord) -> bool:
    """if record has only "N" in seq, return False.

    Args:
        record (SeqRecord.SeqRecord): record

    Returns:
        bool:
    """
    return bool(len(re.sub("[^ATGC]", "", str(record.seq).upper())))


def is_mongrel(record: SeqRecord.SeqRecord) -> bool:
    """Is the creature written in record mongrel ?

    Args:
        record (SeqRecord.SeqRecord): record

    Returns:
        bool:
    """
    name = get_creature_name(record)
    if name is None:
        return False
    else:
        return " x " in name


def is_complete_genome(definition: str) -> bool:
    """完全なミトコンドリアゲノムかどうかを返す.

    Args:
        definition (str): definition

    Returns:
        bool:
    """
    return ("complete" in definition and "shotgun" not in definition
            and "chromosome" not in definition)


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

    if hasattr(record, "annotations"):
        return "contig" in record.annotations.keys()
    else:
        return False


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
    """FileNotFoundError.
    """

    pass
