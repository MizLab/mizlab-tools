# mizlab tools

[![setup.py test](https://github.com/MizLab/mizlab-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/MizLab/mizlab-tools/actions/workflows/ci.yml)


## Installation

```sh
pip install git+https://github.com/MizLab/mizlab-tools.git@main
```

### Uninstall

```sh
pip uninstall mizlab-tools
```

### Usage

#### In `Python` script

```python
import mizlab_tools
```

See in [This docment](https://mizlab.github.io/mizlab-tools/).


#### CLI commands

This library provide some CLI commands. 

- `fetch_gbk`
    ```
    usage: fetch_gbk.py [-h] --email EMAIL [--stdin] [-d DESTINATION]
                        [accessions ...]

    Fetch some gbk files.

    positional arguments:
      accessions            accession numbers

    optional arguments:
      -h, --help            show this help message and exit
      --email EMAIL         your e-mail address to fetch data. (required)
      --stdin               if you want to use stdin, use this flag.
      -d DESTINATION, --destination DESTINATION
                            If set this property, information is stored into
                            ./<-d>, else ./fetched_data/
    ```

    If you exec bellow.

    ```sh
    fetch_gbk  --email YOUR_EMAIL_ADDRESS NC_012920
    ```

    In `./fetched_data`, there will be `NC_012920.gbk`.
- `fetch_taxon`
    ```
    usage: fetch_taxon.py [-h] {ncbi,gnr} ...

    Fetch taxonomy information.

    positional arguments:
      {ncbi,gnr}
        ncbi      Fetch from NCBI.
        gnr       Fetch from Global Names Resolver

    optional arguments:
      -h, --help  show this help message and exit
    ```

    Include 2 subcommands.

    - `ncbi`
        ```
        usage: fetch_taxon.py ncbi [-h] --email EMAIL [--stdin] [-d DESTINATION]
                                   [taxonomy_ids ...]

        positional arguments:
          taxonomy_ids          taxonomy ids

        optional arguments:
          -h, --help            show this help message and exit
          --email EMAIL         your e-mail address to fetch data. (required)
          --stdin               if you want to use stdin, use this flag.
          -d DESTINATION, --destination DESTINATION
                                If set this property, information is stored into
                                ./<-d>, else stdout.
        ```

    - `gnr`
        ```
        usage: fetch_taxon.py gnr [-h] [--stdin] [-d DESTINATION] [binomial_names ...]

        positional arguments:
          binomial_names        binomial names

        optional arguments:
          -h, --help            show this help message and exit
          --stdin               if you want to use stdin, use this flag.
          -d DESTINATION, --destination DESTINATION
                                If set this property, information is stored into
                                ./<-d>, else stdout.
        ```
- `calculate_weights`
    ```
    usage: calculate_weights.py [-h] [-d DESTINATION] [--allow_chars ALLOW_CHARS]
                                gbkfiles [gbkfiles ...]

    Calculate 3 words weight, based on frequency of appearance.

    positional arguments:
      gbkfiles              Genbank format file paths.

    optional arguments:
      -h, --help            show this help message and exit
      -d DESTINATION, --destination DESTINATION
                            if set this option, the weight save into
                            <destination>/weight.json, else $PWD/weight.json.
      --allow_chars ALLOW_CHARS
                            allowed chars.
    ```

- `calculate_coordinates`
    ```
    usage: calculate_coordinates.py [-h] [--weight WEIGHT] --mapping MAPPING
                                    [-d DESTINATION]
                                    gbkfiles [gbkfiles ...]

    Calculate dna data to coordinate. output stdout or file

    positional arguments:
      gbkfiles              GBK format file path.

    optional arguments:
      -h, --help            show this help message and exit
      --weight WEIGHT       The json format weight file.
      --mapping MAPPING     The json format vector mapping file.
      -d DESTINATION, --destination DESTINATION
                            if set this option, the data is store to
                            <destination>/<accession>.dat
    ```
