#!/usr/bin/env python3

import argparse
import collections
import json
import os
import re
from typing import Counter, Dict, Iterable, Union

import numpy as np
from Bio import Seq, SeqIO
from nptyping import Float64

from . import gbk_utils


def calc_weights(sequences: Iterable[Union[str, Seq.Seq]],
                 allowed: str) -> Dict[str, Float64]:
    """calculate the weights .

    Args:
        sequences (Iterable[Union[str, Seq.Seq]]): Itatable object like sequence.
        allowed (str): allowed string like "ATGC"

    Returns:
        Dict[str, Float64]:
    """
    counter: Counter[str] = collections.Counter()
    for sequence in sequences:
        counter += count_words(sequence, allowed, 3)

    weights: Dict[str, Float64] = compute_self_entropy(counter, allowed)
    return weights


def count_words(string: Union[str, Seq.Seq], allowed: str, length: int) -> Counter[str]:
    """count words split per {length}

    Args:
        string (Union[str, Seq.Seq]): string
        allowed (str): allowed string like "ATGC"
        length (int): separation length.

    Returns:
        Counter[str]:
    """
    pretty = re.sub(f"[^{allowed}]", "", str(string))
    return Counter(
        map(lambda x: "".join(x), gbk_utils.window_search(pretty, length, None)))


def compute_self_entropy(counter: Counter[str], allowed: str) -> Dict[str, Float64]:
    """compute_self_entropy.

    Args:
        counter (Counter[str]): counter of words.
        allowed (str): allowed string like "ATGC"

    Returns:
        Dict[str, Float64]:
    """
    bases = set(allowed)
    entropies: Dict[str, float] = {}
    for body, content in counter.items():
        head = body[:-1]
        denominator = sum([counter.get(head + base, 0) for base in bases])
        entropies[body] = -1 * np.log2(content / denominator)
    return entropies


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate 3 words weight, based on frequency of appearance.")
    parser.add_argument("gbkfiles", nargs="+", help="Genbank format file paths.")
    parser.add_argument(
        "-d",
        "--destination",
        help=
        "if set this option, the weight save into <destination>/weight.json, else $PWD/weight.json."
    )
    parser.add_argument("--allow_chars", default="ATGC", help="allowed chars.")

    args = parser.parse_args()

    sequences = []
    for p in args.gbkfiles:
        for record in SeqIO.parse(p, "genbank"):
            sequences.append(record.seq)

    weights = calc_weights(sequences, args.allow_chars)
    weight_for_output: Dict[str, float] = {k: float(v) for k, v in weights.items()}

    if args.destination is not None:
        os.makedirs(args.destination, exist_ok=True)
        with open(os.path.join(args.destination, "weights.json"), "w") as f:
            json.dump(weight_for_output, f)
    else:
        print(json.dumps(weight_for_output))


if __name__ == "__main__":
    main()
