#!/usr/bin/env python3

import argparse
import os
import sys
from typing import Iterable, Iterator

from Bio import Entrez, SeqIO, SeqRecord


def fetch(accession_numbers: Iterable[str],
          email: str) -> Iterator[SeqRecord.SeqRecord]:
    Entrez.email = email
    n_once = 1000
    records = tuple(accession_numbers)
    for i in range(0, len(records), n_once):
        with Entrez.efetch(db="nucleotide",
                           id=records[i:i + n_once],
                           rettype="gb",
                           retmode="text") as handle:
            for r in SeqIO.parse(handle, "genbank"):
                yield r


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch some gbk files.")
    parser.add_argument("--email",
                        required=True,
                        help="your e-mail address to fetch data. (required)")
    parser.add_argument("accessions", nargs="*", help="accession numbers")
    parser.add_argument("--stdin",
                        action="store_true",
                        help="if you want to use stdin, use this flag.")
    args = parser.parse_args()

    os.makedirs("fetched_data", exist_ok=True)

    if args.stdin:
        if sys.stdin.isatty:
            n_records = int(input())
            accessions = [input() for _ in range(n_records)]
        else:
            accessions = [line for line in sys.stdin.readlines()]
        for fetched in fetch(accessions, args.email):
            acc = fetched.name
            with open(f"fetched_data/{acc}.gbk", "w") as f:
                SeqIO.write(fetched, f, "genbank")
    else:
        for fetched in fetch(args.accessions, args.email):
            acc = fetched.name
            with open(f"fetched_data/{acc}.gbk", "w") as f:
                SeqIO.write(fetched, f, "genbank")
