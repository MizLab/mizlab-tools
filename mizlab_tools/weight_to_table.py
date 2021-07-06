import argparse
import json
from typing import Dict, List, Union


def weights_to_table(weights: Dict[str, float], firsts: str, seconds: str,
                     thirds: str) -> List[List[Union[str, float]]]:

    # make empty table
    table = [[""] * (len(seconds) + 2) for _ in range(len(firsts) * len(thirds) + 1)]

    for i, s in enumerate(firsts):
        for j, ss in enumerate(thirds):
            table[(i * len(thirds) + 1) + j][0] = s
            table[(i * len(thirds) + 1) + j][-1] = ss
    for i, s in enumerate(seconds):
        table[0][i + 1] = s

    for row in range(1, len(table)):
        for col in range(1, len(seconds) + 1):
            triptlet = table[row][0] + table[0][col] + table[row][-1]
            table[row][col] = weights[triptlet]

    return table


def main():
    parser = argparse.ArgumentParser(description="Convert weights to table like csv.")
    parser.add_argument("source_json", nargs=1, help="The path to source weight json.")
    parser.add_argument("--stdin", action="store_true", help="Use stdin.")
    parser.add_argument("--delimiter",
                        default=",",
                        help="Delimiter for table. Default to ','")
    args = parser.parse_args()

    with open(args.source_json[0]) as f:
        d = json.load(f)

    firsts = set([])
    seconds = set([])
    thirds = set([])

    for key in d.keys():
        firsts.add(key[0])
        seconds.add(key[1])
        thirds.add(key[2])

    table = weights_to_table(
        d,
        "".join(sorted(firsts)),
        "".join(sorted(seconds)),
        "".join(sorted(thirds)),
    )

    for row in table:
        print(args.delimiter.join(map(str, row)))


if __name__ == "__main__":
    main()
