#!/usr/bin/env python3

import argparse
import json
import os
import sys
from collections import deque
from pathlib import Path
from typing import (Any, Deque, Dict, Iterable, Iterator, List, Literal,
                    Optional, Sequence, Tuple, Union)

import requests
from Bio import Entrez, SeqIO, SeqRecord

from . import gbk_utils, utils


def fetch_taxon(records: Iterable[SeqRecord.SeqRecord],
                email: str,
                n_once: int = 100) -> Dict[str, Union[str, Dict[str, Dict[str, str]]]]:
    """Fetch taxonomy information from NCBI and some DB via Global Names Resolver.

    Args:
        records (Iterable[SeqRecord.SeqRecord]): records
        email (str): email address used in fetching Entrez.
        n_once (int): number of taxonomy informations at once fetching.

    Returns:
        Dict[str, Union[str, Dict[str, Dict[str, str]]]]:
    """
    taxon_info = {}
    for s in utils.split_per_n(records, n=n_once):
        subrecords = tuple(s)
        taxon_ids = filter(lambda x: x is not None,
                           map(gbk_utils.get_taxonID, subrecords))
        binomial_names = tuple(
            filter(lambda x: x is not None, map(gbk_utils.get_creature_name,
                                                subrecords)))
        sub_results = {}
        for accession, binomial_name, res_ncbi, res_gnr in zip(
                map(lambda x: x.name, subrecords), binomial_names,
                fetch_taxon_from_NCBI(taxon_ids, email=email, n_once=n_once),
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
        n_once: int = 100,
        verbose: bool = False) -> Iterator[Union[List[Dict[str, str]], Dict[str, str]]]:
    """Fetch taxonomy informations from NCBI(Entrez).

    Args:
        taxonomy_ids (Iterable[str]): taxonomy_ids
        email (str): email address used in fetching Entrez.
        n_once (int): number of taxonomy informations at once fetching.
        verbose (bool): whether include ambigous taxonomy informations.

    Returns:
        Iterator[Union[List[Dict[str, str]], Dict[str, str]]]:
    """
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


def fetch_taxon_from_GNR(names: Iterable[str],
                         priority: Optional[Sequence[int]] = None,
                         n_once: int = 100) -> Iterator[Dict[str, Dict[str, str]]]:
    """fetch_taxon_from_GNR.

    Args:
        names (Iterable[str]): names in biology.
        priority (Optional[Sequence[int]]): priority. see in https://resolver.globalnames.org
        n_once (int): number of taxonomy informations at once fetching.

    Returns:
        Iterator[Dict[str, Dict[str, str]]]:
    """
    url = "http://resolver.globalnames.org/name_resolvers.json"
    if priority is None:
        priority = [4, 3, 179, 11, 1, 8]
    for s in utils.split_per_n(names, n_once):
        queue: Deque = deque()
        queue.append(tuple(s))
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
            else:    # when not Exception
                for entry in data["data"]:
                    yield _parse_gnr_response(entry)


def _parse_gnr_response(entry: Dict[str, Any]) -> Dict[str, Dict[str, str]]:
    """Parse responce from Grobal Names Resolver.

    Args:
        entry (Dict[str, Any]): The responce of a query from global names resolver.

    Returns:
        Dict[str, Dict[str, str]]:
    """
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
    """Wrapper of fetch_taxon_from_NCBI.

    Args:
        args (argparse.Namespace): arguments from commandline.

    Returns:
        None:
    """
    dst = Path()
    if args.destination is not None:
        dst = Path(args.description)
        dst.mkdir(parents=True, exist_ok=True)

    if args.gbk:
        output_names, taxon_ids = _filename_and_search_key(args.taxonomy_ids,
                                                           "taxon_id")
    else:
        output_names, taxon_ids = zip(args.taxon_ids, args.taxon_ids)

    for taxonid, record in zip(args.taxonomy_ids,
                               fetch_taxon_from_NCBI(args.taxonomy_ids, args.email)):
        if args.destination is not None:
            with open(str(dst / f"{taxonid}.json"), "w") as f:
                json.dump(record, f)
        else:
            print(json.dumps(record, indent=2))


def _fetch_gnr(args: argparse.Namespace) -> None:
    """Wrapper of fetch_taxon_from_GNR

    Args:
        args (argparse.Namespace): arguments from commandline.

    Returns:
        None:
    """
    dst = Path()
    if args.destination is not None:
        dst = Path(args.destination)
        dst.mkdir(parents=True, exist_ok=True)

    if args.gbk:
        output_names, binomial_names = _filename_and_search_key(
            args.binomial_names, "binomial_name")
    else:
        output_names, binomial_names = zip(args.binomial_names, args.binomial_names)

    for output_name, record in zip(output_names, fetch_taxon_from_GNR(binomial_names)):
        if args.destination is not None:
            with open(str(dst / f"{output_name}.json"), "w") as f:
                json.dump(record, f)
        else:
            print(json.dumps(record, indent=2))


def _filename_and_search_key(
    gbk: gbk_utils.Openable,
    search_key: Literal["taxon_id", "binomial_name"],
) -> Iterator[Tuple[str, str]]:
    for record in SeqIO.parse(gbk, "genbank"):
        if search_key == "taxon_id":
            key = gbk_utils.get_taxonID(record)
            if key is None:
                raise ValueError
            yield record.name, key
        else:
            key = gbk_utils.get_creature_name(record)
            if key is None:
                raise ValueError
            yield record.name, key


def _fetch_all(args: argparse.Namespace) -> None:
    records = []
    for f in args.gbkfiles:
        for record in SeqIO.parse(f, "genbank"):
            records.append(record)
    results = fetch_taxon(records, email=args.email)

    dst = Path()
    if args.destination is not None:
        dst = Path(args.destination)
        dst.mkdir(parents=True, exist_ok=True)
        with open(dst / "result.json", "w") as f:
            json.dump(results, f)
    else:
        print(json.dumps(results, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch taxonomy information.")
    subparser = parser.add_subparsers()
    parser_ncbi = subparser.add_parser("ncbi", help="Fetch from NCBI.")
    parser_ncbi.add_argument("--email",
                             required=True,
                             help="your e-mail address to fetch data. (required)")
    parser_ncbi.add_argument("--stdin", action="store_true", help="use stdin.")
    parser_ncbi.add_argument("--gbk",
                             action="store_true",
                             help="use Genbank form data.")
    parser_ncbi.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser_ncbi.add_argument("taxonomy_ids", nargs="+", help="taxonomy ids")
    parser_ncbi.set_defaults(hander=_fetch_ncbi)

    parser_gnr = subparser.add_parser("gnr", help="Fetch from Global Names Resolver")
    parser_gnr.add_argument("--stdin", action="store_true", help="use stdin.")
    parser_gnr.add_argument("--gbk", action="store_true", help="use Genbank form data.")
    parser_gnr.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser_gnr.add_argument("binomial_names", nargs="+", help="binomial names")
    parser_gnr.set_defaults(hander=_fetch_gnr)

    parser_all = subparser.add_parser("all", help="Fetch from NCBI and GNR.")
    parser_all.add_argument("--email",
                            required=True,
                            help="your e-mail address to fetch data. (required)")
    parser_all.add_argument("--stdin", action="store_true", help="use stdin.")
    parser_all.add_argument(
        "-d",
        "--destination",
        help="If set this property, information is stored into ./<-d>, else stdout.")
    parser_all.add_argument("gbkfiles", nargs="+", help="Path of Genbank from data.")
    parser_all.set_defaults(hander=_fetch_all)

    args = parser.parse_args()
    if hasattr(args, "hander"):
        args.hander(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
