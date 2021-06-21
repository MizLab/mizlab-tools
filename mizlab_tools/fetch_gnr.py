#!/usr/bin/env python3

import argparse
import json
import os
import sys
from typing import Dict, Iterable, Iterator, List, Optional, Tuple, Union

from Bio import Entrez


def fetch_taxons(
        taxonomy_ids: Union[List[str], Tuple[str]],
        email: str,
        verbose: bool = False) -> Iterator[Union[List[Dict[str, str]], Dict[str, str]]]:
    Entrez.email = email
    n_once = 100

    for i in range(0, len(taxonomy_ids), n_once):
        with Entrez.efetch(db="Taxonomy", id=taxonomy_ids[i:i + n_once]) as handle:
            records = Entrez.read(handle)
            for record in records:
                if verbose:
                    yield [{
                        "ScientificName": lineage["ScientificName"],
                        "Rank": lineage["Rank"].lower()
                    } for lineage in record["LineageEx"]]
                else:
                    yield {
                        lineage["Rank"].lower(): lineage["ScientificName"]
                        for lineage in record["LineageEx"]
                        if lineage["Rank"] != "no rank"
                    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch taxonomy information from Entrez.")
    parser.add_argument("--email",
                        required=True,
                        help="your e-mail address to fetch data. (required)")
    parser.add_argument("--stdin",
                        action="store_true",
                        help="if you want to use stdin, use this flag.")
    parser.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser.add_argument("--verbose",
                        action="store_true",
                        help="Show ambiguous taxon info.")
    parser.add_argument("taxonomy_ids", nargs="*", help="taxonomy ids")
    args = parser.parse_args()

    if args.destination:
        os.makedirs(args.destination)

    if args.stdin:
        if sys.stdin.isatty():
            n_inputs = int(input())
            entries = [input() for _ in range(n_inputs)]
        else:
            entries = [line for line in sys.stdin.readlines()]
    else:
        entries = args.taxonomy_ids

    for taxon_id, result in zip(entries, fetch_taxon(entries, args.email,
                                                     args.verbose)):
        if not args.destination:
            print(f"---{taxon_id}---")
            print(result)
        else:
            with open(f"{args.destination}/{taxon_id}.json", "w") as f:
                json.dump(result, f, indent=2)
