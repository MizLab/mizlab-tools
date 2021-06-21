#!/usr/bin/env python3

import argparse
import json
import os
import re
from collections import Counter
from itertools import product
from typing import Any, Dict, List, Union

from Bio import SeqIO

Openable = Union[str, bytes, int, "os.PathLike[Any]"]


def calculate_weights(paths: List[Openable], allow_chars: str) -> Dict[str, float]:
    triplet_counter: Counter[str] = Counter()
    for p in paths:
        for record in SeqIO.parse(p, "genbank"):
            seq = re.sub(f"[^{allow_chars}]", "", str(record.seq))
            c = Counter([seq[i:i + 3] for i in range(len(seq) - 2)])
            triplet_counter += c

    weights = {}
    for first, second in product(allow_chars, repeat=2):
        twins = first + second
        sum_of_start_with_twins = sum([triplet_counter[twins + t] for t in allow_chars])
        for third in allow_chars:
            codon = twins + third
            weights[codon] = triplet_counter[codon] / sum_of_start_with_twins
    return weights


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate 3 words weight, based on frequency of appearance.")
    parser.add_argument("gbkfiles", nargs="+", help="Genbank format file paths.")
    parser.add_argument(
        "--stdout",
        action="store_true",
        help="if set this flag, show weight into stdout as json-format.")
    parser.add_argument(
        "-d",
        "--destination",
        help=
        "if set this option, the weight save into <destination>/weight.json, else $PWD/weight.json."
    )
    parser.add_argument("--allow_chars", default="ATGC", help="allowed chars.")

    args = parser.parse_args()

    weights = calculate_weights(args.gbkfiles, allow_chars=args.allow_chars)
    if args.destination is not None:
        os.makedirs(args.destination, exist_ok=True)
        with open(os.path.join(args.destination, "weights.json"), "w") as f:
            json.dump(weights, f)
    else:
        print(json.dumps(weights))
