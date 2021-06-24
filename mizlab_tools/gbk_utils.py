import os
import re
from collections import Counter
from pathlib import Path
from typing import (Any, AnyStr, Dict, Iterable, Iterator, Optional, Sequence, Tuple,
                    TypeVar, Union)

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


def has_contig(record: SeqRecord.SeqRecord) -> bool:
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
