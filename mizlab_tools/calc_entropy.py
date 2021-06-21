import re
from collections import Counter
from typing import AnyStr, Dict, Iterable

import numpy as np

from . import gbk_utils


def calc_weights(sequences: Iterable[AnyStr],
                 atgc_only: bool = False) -> Dict[str, float]:
    """塩基配列のリストから3塩基ごとの情報量を算出し、返す
       返す配列の塩基文字は大文字である

    Args:
        sequences (Iterable[AnyStr]): 配列の集合
        atgc_only (bool): Defaults to False.

    Returns:
        Dict[str, float]: 3塩基に対する情報量の辞書
    """
    counter = count_triplets(map(lambda x: x.upper(), sequences), atgc_only=atgc_only)
    return compute_self_entropy(counter)


def count_triplets(sequences: Iterable[AnyStr], atgc_only: bool = False) -> Counter:
    """塩基配列の集合から3文字の出現個数をカウントし、返す

    Args:
        sequences (Iterable[AnyStr]): 配列の集合
        atgc_only (bool): Defaults to False.

    Returns:
        Counter: 3塩基に対する出現数の辞書
    """
    counter = Counter()
    for seq in sequences:
        if atgc_only:
            seq = re.sub("[^ATGC]", "", seq)
        c = Counter(
            map(lambda x: "".join(x), gbk_utils.window_search(seq, window_size=3)))
        counter += c
    return counter


def compute_self_entropy(triplets: Counter) -> dict:
    """自己情報量に基づいて3塩基のdictのvalueをエントロピーに書き換える

    Args:
        triplets (Counter): 3塩基に対する出現数の辞書

    Returns:
        dict: 3塩基に対する情報量の辞書
    """
    bases = {"A", "T", "G", "C"}
    entropies = {}
    for triplet, content in triplets.items():
        first_2_base = triplet[:2]
        denominator = sum([triplets[first_2_base + base] for base in bases])
        entropies[triplet] = -1 * np.log2(content / denominator)
    return entropies
