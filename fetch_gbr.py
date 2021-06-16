#!/usr/bin/env python3

import argparse
import json
import os
import sys
from typing import (Any, AnyStr, Dict, Iterable, Iterator, List, Optional,
                    Tuple, Union)

import requests


def fetch_taxon(
    names: Union[List[str], Tuple[str]],
    priority: Optional[Iterable[Union[str,
                                      int]]]) -> Iterator[Dict[str, Dict[str, str]]]:
    invalid_ranks_and_names = {"", "no rank", "no rank - terminal", "Not assingned"}
    if priority is None:
        priority = [4, 3, 179, 11, 1, 8]    # decide by heuristic method
    results: Dict[str, Dict[str, str]] = {}
    base_params: Dict[str, Union[bool, str]] = {
        "best_match_only": True,
        "preferred_data_sources": "|".join(map(str, priority))
    }
    url = "http://resolver.globalnames.org/name_resolvers.json"

    for i in range(0, len(names), 100):
        p: Dict[str, Union[bool, str]] = {"names": "|".join(names[i:i + 100])}
        p.update(base_params)
        r = requests.get(url, params=p)
        data = r.json()
        for entry in data["data"]:
            results = {}
            for record in entry.get("preferred_results", []):
                parts = {}
                for rank, name in zip(record["classification_path_ranks"].split("|"),
                                      record["classification_path"].split("|")):
                    if (rank in invalid_ranks_and_names
                            or name in invalid_ranks_and_names):
                        continue
                    else:
                        parts[rank] = name
                results[record["data_source_title"]] = parts
            yield results


def _parse_gnr_response(entry: Dict[str, Any]) -> Dict[str, Union[str, Dict[str, str]]]:
    """_parse_gnr_response.

    Args:
        entry (GnrResponse): entry

    Returns:
        Dict[str, Union[str, Dict[str, str]]]: Like bellow.
            "homo sapiens": {
                "NCBI": {"Class": "Mammalia"...},
                ...
             }
    """
    invalid_ranks_and_names = {"", "no rank", "no rank - terminal", "Not assingned"}
    results: Dict[str, Union[str, Dict[str, str]]] = {}
    for record in entry.get("preferred_results", []):
        ind_db = {}
        for rank, name in zip(record["classification_path_ranks"].split("|"),
                              record["classification_path"].split("|")):
            if rank in invalid_ranks_and_names or name in invalid_ranks_and_names:
                continue
            else:
                ind_db[rank] = name
        results[record["data_source_title"]] = ind_db
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch taxonomy information from Grobal Names Resolver.")
    parser.add_argument("--stdin",
                        action="store_true",
                        help="if you want to use stdin, use this flag.")
    parser.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser.add_argument("names", nargs="*", help="names in biology.")
    parser.add_argument("-p", "--priority", nargs="*", help="Priority of DB.")
    args = parser.parse_args()

    if args.destination:
        os.makedirs(args.destination, exist_ok=True)

    if args.stdin:
        if sys.stdin.isatty():
            n_inputs = int(input())
            entries = [input() for _ in range(n_inputs)]
        else:
            entries = [line for line in sys.stdin.readlines()]
    else:
        entries = args.names

    for name, record in zip(entries, fetch_taxon(entries, args.priority)):
        if not args.destination:
            print(f"---{name}---")
            print(record)
        else:
            with open(f"{args.destination}/{name}.json", "w") as f:
                json.dump(record, f, indent=2)
