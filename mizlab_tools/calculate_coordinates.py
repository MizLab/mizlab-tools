#!/usr/bin/env python3

import argparse
import json
import os
import re
from typing import Dict, List, Optional, Union

from Bio import SeqIO, SeqRecord

from . import gbk_utils


def calc_coord(
        record: SeqRecord.SeqRecord,
        mapping: Dict[str, List[int]],
        weight: Optional[Dict[str, Union[int, float]]] = None) -> List[List[float]]:
    """calculate coordinates with mapping of vencor to base.

    Args:
        record (SeqRecord.SeqRecord): record
        mapping (Dict[str, List[int]]): mapping of vector
        weight (Optional[Dict[str, Union[int, float]]]): weight of some base combi.

    Returns:
        List[List[float]]: result of calculation.
    """
    allow_str = "".join(mapping.keys())

    filterd_seq = re.sub(f"[^{allow_str}]", "", str(record.seq))

    dimension = 0    # default value, it will be rewrited
    for v in mapping.values():
        dimension = len(v)
        break

    if weight is None:
        weight = {}

    coordinates = [[0. for _ in range(dimension)]]
    for window in gbk_utils.window_search(filterd_seq, 3, overhang="before"):
        triplet = "".join(map(str, window))
        rate = weight.get(triplet, 1.0)
        vector = list(map(lambda x: x * rate, mapping[triplet[-1]]))
        coordinates.append([a + b for a, b in zip(coordinates[-1], vector)])

    return coordinates


def main() -> None:
    """main.

    Args:

    Returns:
        None:
    """
    parser = argparse.ArgumentParser(
        description="Calculate dna data to coordinate. output stdout or file")
    parser.add_argument("gbkfiles", nargs="+", help="GBK format file path.")
    parser.add_argument("--weight", help="The json format weight file.")
    parser.add_argument("--mapping",
                        help="The json format vector mapping file.",
                        required=True)
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
            data = calc_coord(record, mapping, weight)
            if destination is not None:
                with open(f"{destination}/{record.name}.dat", "w") as f:
                    for d in data:
                        print(" ".join(map(str, d)), file=f)
            else:
                for d in data:
                    print(" ".join(map(str, d)))


if __name__ == "__main__":
    main()
