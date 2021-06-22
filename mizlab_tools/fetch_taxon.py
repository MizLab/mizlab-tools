#!/usr/bin/env python3

import argparse
import json
import os
import sys
from collections import deque
from pathlib import Path
from typing import (Any, Deque, Dict, Iterable, Iterator, List, Optional,
                    Sequence, Tuple, Union)
from urllib.parse import ParseResult

import requests
from Bio import Entrez, SeqRecord

from . import gbk_utils, utils


def fetch_taxon(records: Iterable[SeqRecord.SeqRecord],
                email: str,
                n_once: int = 100) -> Dict[str, Union[str, Dict[str, Dict[str, str]]]]:

    taxon_info = {}
    for s in utils.split_per_n(records, n=n_once):
        subrecords = tuple(records)
        taxon_ids = filter(lambda x: x is not None,
                           map(gbk_utils.get_taxonID, subrecords))
        binomial_names = tuple(
            filter(lambda x: x is not None, map(gbk_utils.get_creature_name,
                                                subrecords)))
        sub_results = {}
        for accession, binomial_name, res_ncbi, res_gnr in zip(
                map(lambda x: x.name, subrecords), binomial_names,
                fetch_taxon_from_NCBI(taxon_ids, email=email),
                fetch_taxon_from_GNR(binomial_names)):
            sub_results[accession] = {
                "binomial_name": binomial_name,
                "taxon": {
                    **{
                        "NCBI Taxonomy": res_ncbi.copy()
                    }
                },
                **res_gnr
            }
        taxon_info.update(sub_results)
    return taxon_info


def fetch_taxon_from_NCBI(
        taxonomy_ids: Iterable[str],
        email: str,
        verbose: bool = False) -> Iterator[Union[List[Dict[str, str]], Dict[str, str]]]:
    Entrez.email = email
    n_once = 100

    for subset in utils.split_per_n(taxonomy_ids, n=n_once):
        with Entrez.efetch(db="Taxonomy", id=tuple(subset)) as handle:
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


def fetch_taxon_from_GNR(
        names: Iterable[str],
        priority: Optional[Sequence[int]] = None
) -> Iterator[Dict[str, Dict[str, str]]]:
    url = "http://resolver.globalnames.org/name_resolvers.json"
    if priority is None:
        priority = [4, 3, 179, 11, 1, 8]
    queue: Deque = deque()
    queue.append(tuple(names))
    while (queue):
        sub_names = queue.popleft()
        params = {
            "names": "|".join(sub_names),
            "best_match_only": True,
            "preferred_data_sources": "|".join(map(str, priority))
        }
        r = requests.get(url, params=params)
        try:
            data = r.json()
        except Exception:
            if len(sub_names) == 1:
                yield {}
            else:
                mid = len(sub_names) // 2
                queue.append(sub_names[:mid])
                queue.append(sub_names[mid:])
        else:
            for entry in data["data"]:
                yield _parse_gnr_response(entry)


def _parse_gnr_response(entry: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    invalid_ranks_and_names = {"", "no rank", "no rank - terminal", "Not assingned"}
    results = {}
    for record in entry.get("preferred_results", []):
        db_res = {}
        for rank, name in zip(record["classification_path_ranks"].split("|"),
                              record["classification_path"].split("|")):
            if rank in invalid_ranks_and_names or name in invalid_ranks_and_names:
                continue
            else:
                db_res[rank] = name
        results[record["data_source_title"]] = db_res
    return results


def _fetch_ncbi(args: argparse.Namespace) -> None:
    dst = Path()
    if args.destination is not None:
        dst = Path(args.description)
        dst.mkdir(parents=True, exist_ok=True)

    for taxonid, record in zip(args.taxonomy_ids,
                               fetch_taxon_from_NCBI(args.taxonomy_ids, args.email)):
        if args.destination is not None:
            with open(str(dst / f"{taxonid}.json"), "w") as f:
                json.dump(record, f)
        else:
            print(json.dumps(record, indent=2))


def _fetch_gnr(args: argparse.Namespace) -> None:
    dst = Path()
    if args.destination is not None:
        dst = Path(args.destination)
        dst.mkdir(parents=True, exist_ok=True)
    for binomial_name, record in zip(args.binomial_names,
                                     fetch_taxon_from_GNR(args.binomial_names)):
        if args.destination is not None:
            with open(str(dst / f"{binomial_name}.json"), "w") as f:
                json.dump(record, f)
        else:
            print(json.dumps(record, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch taxonomy information.")
    subparser = parser.add_subparsers()
    parser_ncbi = subparser.add_parser("ncbi", help="Fetch from NCBI.")
    parser_ncbi.add_argument("--email",
                             required=True,
                             help="your e-mail address to fetch data. (required)")
    parser_ncbi.add_argument("--stdin",
                             action="store_true",
                             help="if you want to use stdin, use this flag.")
    parser_ncbi.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser_ncbi.add_argument("taxonomy_ids", nargs="*", help="taxonomy ids")
    parser_ncbi.set_defaults(hander=_fetch_ncbi)

    parser_gnr = subparser.add_parser("gnr", help="Fetch from Global Names Resolver")
    parser_gnr.add_argument("--stdin",
                            action="store_true",
                            help="if you want to use stdin, use this flag.")
    parser_gnr.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser_gnr.add_argument("binomial_names", nargs="*", help="binomial names")
    parser_gnr.set_defaults(hander=_fetch_gnr)

    args = parser.parse_args()
    if hasattr(args, "hander"):
        args.hander(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
    # parser = argparse.ArgumentParser(
    #     description="Fetch taxonomy information from Entrez.")
    # parser.add_argument("--email",
    #                     required=True,
    #                     help="your e-mail address to fetch data. (required)")
    # parser.add_argument("--stdin",
    #                     action="store_true",
    #                     help="if you want to use stdin, use this flag.")
    # parser.add_argument(
    #     "-d",
    #     "--destination",
    #     help="If set this property, information is stored into ./<-d>, else stdout.")
    # parser.add_argument("--verbose",
    #                     action="store_true",
    #                     help="Show ambiguous taxon info.")
    # parser.add_argument("taxonomy_ids", nargs="*", help="taxonomy ids")
    # args = parser.parse_args()
    #
    # if args.destination:
    #     os.makedirs(args.destination)
    #
    # if args.stdin:
    #     if sys.stdin.isatty():
    #         n_inputs = int(input())
    #         entries = [input() for _ in range(n_inputs)]
    #     else:
    #         entries = [line for line in sys.stdin.readlines()]
    # else:
    #     entries = args.taxonomy_ids
    #
    # for taxon_id, result in zip(entries, fetch_taxon(entries, args.email,
    #                                                  args.verbose)):
    #     if not args.destination:
    #         print(f"---{taxon_id}---")
    #         print(result)
    #     else:
    #         with open(f"{args.destination}/{taxon_id}.json", "w") as f:
    #             json.dump(result, f, indent=2)
