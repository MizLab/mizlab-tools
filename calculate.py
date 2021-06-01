#!/usr/bin/env python3

import argparse
import json
import os
import re
from typing import (Dict, Iterable, Iterator, List, Literal, Optional, Tuple,
                    TypeVar, Union)

from Bio import SeqIO, SeqRecord

OVERHANG = Literal["before", "after", "both"]
T = TypeVar("T")


def window_search(target: Iterable[T],
                  size: int,
                  overhang: Optional[OVERHANG] = None) -> Iterator[Tuple[T, ...]]:
    fixed_target = tuple(target)

    if overhang == "before" or overhang == "both":
        for i in range(1, size):
            yield fixed_target[:i]

    for i in range(len(fixed_target) - size + 1, ):
        yield fixed_target[i:i + size]

    if overhang == "after" or overhang == "both":
        for i in range(len(fixed_target) - size + 1, len(fixed_target)):
            yield fixed_target[i:]


def factory(record: SeqRecord.SeqRecord,
            mapping: Dict[str, List[int]],
            weight: Optional[Dict[str, Union[int, float]]] = None) -> List[List[float]]:
    allow_str = "".join(mapping.keys())

    filterd_seq = re.sub(f"[^{allow_str}]", "", record.seq)

    dimension = 0    # default value, it will be rewrited
    for v in mapping.values():
        dimension = len(v)
        break

    if weight is None:
        weight = {}

    coordinates = [[0. for _ in range(dimension)]]
    for t in window_search(filterd_seq, 3, overhang="both"):
        triplet = "".join(t)
        rate = weight.get(triplet, 1.0)
        vector = list(map(lambda x: x * rate, mapping[triplet[-1]]))
        coordinates.append([a + b for a, b in zip(coordinates[-1], vector)])

    return coordinates


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate dna data to coordinate. output stdout or file")
    parser.add_argument("gbkfiles", nargs="+", help="GBK format file path.")
    parser.add_argument("--weight", help="The json format weight file.")
    parser.add_argument("--mapping", help="The json format vector mapping file.")
    parser.add_argument(
        "-d",
        "--destination",
        help="if set this option, the data is store to <destination>/<accession>.dat")
    args = parser.parse_args()

    if args.weight is not None:
        with open(args.weight) as f:
            weight = json.load(f)
    else:
        weight = None

    if args.destination is None:
        destination = None    # stdout
    else:
        destination = args.destination
        os.makedirs(destination, exist_ok=True)

    with open(args.mapping) as f:
        mapping = json.load(f)

    if args.destination is not None:
        os.makedirs(args.destination, exist_ok=True)
    for gbk in args.gbkfiles:
        for record in SeqIO.parse(gbk, "genbank"):
            data = factory(record, mapping, weight)
            if destination is not None:
                with open(f"{destination}/{record.name}.dat", "w") as f:
                    for d in data:
                        print(" ".join(map(str, d)), file=f)
            else:
                for d in data:
                    print(" ".join(map(str, d)))
