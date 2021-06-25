#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import Iterable, Iterator

from Bio import Entrez, SeqIO, SeqRecord

from . import utils


def fetch(accession_numbers: Iterable[str],
          email: str) -> Iterator[SeqRecord.SeqRecord]:
    """Fetch gbk files.

    Args:
        accession_numbers (Iterable[str]): acsession numbers.
        email (str): email used in fetch NCBI(Entrez).

    Returns:
        Iterator[SeqRecord.SeqRecord]:
    """
    Entrez.email = email
    n_once = 1000
    for accs in utils.split_per_n(accession_numbers, n=n_once):
        with Entrez.efetch(db="nucleotide", id=accs, rettype="gb",
                           retmode="text") as handle:
            for r in SeqIO.parse(handle, "genbank"):
                yield r


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch some gbk files.")
    parser.add_argument("--email",
                        required=True,
                        help="your e-mail address to fetch data. (required)")
    parser.add_argument("accessions", nargs="*", help="accession numbers")
    parser.add_argument("--stdin",
                        action="store_true",
                        help="if you want to use stdin, use this flag.")
    parser.add_argument(
        "-d",
        "--destination",
        help=
        "If set this property, information is stored into ./<-d>, else ./fetched_data/")
    args = parser.parse_args()

    if args.destination is None:
        dst = Path("fetched_data")
    else:
        dst = Path(args.destination)

    if not dst.exists():
        dst.mkdir(parents=True, exist_ok=True)

    if args.stdin:
        if sys.stdin.isatty:
            n_records = int(input())
            accessions = [input() for _ in range(n_records)]
        else:
            accessions = [line for line in sys.stdin.readlines()]
        for fetched in fetch(accessions, args.email):
            acc = fetched.name
            with open(str(dst / f"{acc}.gbk"), "w") as f:
                SeqIO.write(fetched, f, "genbank")
    else:
        for fetched in fetch(args.accessions, args.email):
            acc = fetched.name
            with open(str(dst / f"{acc}.gbk"), "w") as f:
                SeqIO.write(fetched, f, "genbank")


if __name__ == "__main__":
    main()
