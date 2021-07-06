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


See in [This docment](https://mizlab.github.io/mizlab-tools/).


#### CLI commands

This library provide some CLI commands.

- `fetch_gbk`

    Fetch Genbank format data.

    If you exec

    ```sh
    fetch_gbk --email YOUR_EMAIL_ADDRESS NC_012920
    ```

    In `./fetched_data`, there will be `NC_012920.gbk`.

- `fetch_taxon`

    Fetch taxonomy information from NCBI | Global Names Resolver.

    Include 3 subcommands.

    - `ncbi`

        if you exec

        ```sh
        fetch_gbk ncbi --email YOUR_EMAIL_ADRESS 9606
        ```

        show taxon information like

        ```json
        {
          "superkingdom": "Eukaryota",
          "clade": "Boreoeutheria",
          "kingdom": "Metazoa",
          "phylum": "Chordata",
          "subphylum": "Craniata",
          "superclass": "Sarcopterygii",
          "class": "Mammalia",
          "superorder": "Euarchontoglires",
          "order": "Primates",
          "suborder": "Haplorrhini",
          "infraorder": "Simiiformes",
          "parvorder": "Catarrhini",
          "superfamily": "Hominoidea",
          "family": "Hominidae",
          "subfamily": "Homininae",
          "genus": "Homo"
        }
        ```

        the command equal
        ```sh
        echo 9606 | fetch_gbk ncbi --email YOUR_EMAIL_ADDRESS --stdin
        ```

        when add `--destination foo` to command, the json dumped at `./foo/9606.json`.

        if 2 or more data given, the number of dumped file will increase accordingly.

        when add flag `--gbk`, you can use `~.gbk` as arguments.

    - `gnr`

        if you exec

        `fetch_taxon gnr "Homo Sapiens"`

        show taxon information on console.

        ```json
        {
          "National Center for Biotechnology Information": {
            "superkingdom": "Eukaryota",
            "clade": "Boreoeutheria",
            "kingdom": "Metazoa",
            "phylum": "Chordata",
            "subphylum": "Craniata",
            "superclass": "Sarcopterygii",
            "class": "Mammalia",
            "superorder": "Euarchontoglires",
            "order": "Primates",
            "suborder": "Haplorrhini",
            "infraorder": "Simiiformes",
            "parvorder": "Catarrhini",
            "superfamily": "Hominoidea",
            "family": "Hominidae",
            "subfamily": "Homininae",
            "genus": "Homo",
            "species": "Homo sapiens"
          },
          "Integrated Taxonomic Information SystemITIS": {
            "Kingdom": "Animalia",
            "Subkingdom": "Bilateria",
            "Infrakingdom": "Deuterostomia",
            "Phylum": "Chordata",
            "Subphylum": "Vertebrata",
            "Infraphylum": "Gnathostomata",
            "Superclass": "Tetrapoda",
            "Class": "Mammalia",
            "Subclass": "Theria",
            "Infraclass": "Eutheria",
            "Order": "Primates",
            "Suborder": "Haplorrhini",
            "Infraorder": "Simiiformes",
            "Superfamily": "Hominoidea",
            "Family": "Hominidae",
            "Subfamily": "Homininae",
            "Genus": "Homo",
            "Species": "Homo sapiens"
          },
          "Open Tree of Life Reference Taxonomy": {
            "domain": "Eukaryota",
            "kingdom": "Metazoa",
            "phylum": "Chordata",
            "subphylum": "Vertebrata",
            "superclass": "Tetrapoda",
            "class": "Mammalia",
            "subclass": "Theria",
            "superorder": "Euarchontoglires",
            "order": "Primates",
            "suborder": "Haplorrhini",
            "infraorder": "Simiiformes",
            "parvorder": "Catarrhini",
            "superfamily": "Hominoidea",
            "family": "Hominidae",
            "subfamily": "Homininae",
            "genus": "Homo",
            "species": "Homo sapiens"
          },
          "GBIF Backbone Taxonomy": {
            "kingdom": "Animalia",
            "phylum": "Chordata",
            "class": "Mammalia",
            "order": "Primates",
            "family": "Hominidae",
            "genus": "Homo",
            "species": "Homo sapiens"
          },
          "Catalogue of Life - June 2021": {
            "unranked": "Biota",
            "kingdom": "Animalia",
            "phylum": "Chordata",
            "class": "Mammalia",
            "subclass": "Theria",
            "infraclass": "Eutheria",
            "order": "Primates",
            "suborder": "Haplorrhini",
            "infraorder": "Simiiformes",
            "superfamily": "Hominoidea",
            "family": "Hominidae",
            "subfamily": "Homininae",
            "genus": "Homo",
            "species": "Homo sapiens"
          },
          "The Interim Register of Marine and Nonmarine Genera": {
            "kingdom": "Animalia",
            "phylum": "Chordata",
            "class": "Mammalia",
            "order": "Primates",
            "family": "Hominidae",
            "genus": "Homo",
            "species": "Homo sapiens"
          }
        }
        ```

        the command equal

        ```sh
        echo "Homo Sapiens" | fetch_taxon gnr
        ```

        when add `--destination foo`, the json dumped at `./foo/"Homo Sapiens.json"`

        if 2 or more data given, the number of dumped file will increase accordingly.

        when add flag `--gbk`, you can use `~.gbk` as arguments.

    - `all`

        This command is a combination of the above 2.

        ```sh
        fetch_taxon all --email YOUR_EMAIL_ADDRESS NC_012920.gbk
        ```

        show the joined json of each one.

        this command allow only path of genbank file as arguments.

- `calculate_weights`

    Calculate weight of 3 bases.

    ```sh
    calculate_weights path/to/gbk/file ...
    ```

    show calculated weights like

    ```json
    {
        "AAA": 1.5,
        ...
        "CCC": 4.2
    }
    ```

    if `--allow_chars` option is given.

    if you want to use M in addition to ATGC, use `--allowed_chars ATGCR`

    `--destination foo` is given, calculated weights dumped at `./foo/weights.json`.

- `calculate_coordinates`

    Translate DNA to graph coordinates.

    ```sh
    calculate_coordinates --mapping path/to/mapping.json --weight path/to/weights.json path/to/gbkfile ...
    ```

    after exec, show coordinates in stdout.

    if you want to dump coordinates, use `--destination`.

    ```sh
    calculate_coordinates --mapping path/to/mapping.json --weight path/to/weights.json --destination ./outputs NC_012920.gbk
    ```

    the coordinate dumped at `./outputs/NC_012920.dat`




## Requirements
    - `Python` >= 3.6.1
