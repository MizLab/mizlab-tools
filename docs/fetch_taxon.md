# fetch_taxon.py

taxonomy IDを元に分類情報を取得するプログラム

## 使い方

```console
$ python fetch_taxon.py --email taxonomy_ids...
```

- `--email`
    `biopython`で`Entrez`を取得する際に必要になるメールアドレス（大学のものでよい）
- `--stdin`
    標準入力でtaxonomy idを受け取るときに指定するフラグ。
    パイプを使わないで渡す場合は入力個数おｗ先に指定する必要がある。
    例
    ```console  
    # パイプ
    $ cat accession.txt | python fetch_taxon.py --email ~ --stdin

    # 非パイプ
    $ python fetch_taxon.py --email ~ --stdin
    1
    6906
    ```
- `-d`,`--destination`
    取得した分類情報をどこのディレクトリに出力するか。
    これが指定されていない時は分類情報を標準出力に書き出す。
- `--verbose`
    このフラグが指定された時、曖昧な分類情報（`class`などの階層名がない）も取得する。

## 動作条件
- python3.x
- biopython
