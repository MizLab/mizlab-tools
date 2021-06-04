# fetch_gbk.py

## 使い方

```console
$ python fetch_gbk.py --email ~ accessions...
```

- `--email`
    `biopython`で`Entrez`を取得する際に必要になるメールアドレス（大学のものでよい）
- `--stdin`
    標準入力でアクセッション番号を受け取るときに指定するフラグ
    パイプを使わないで渡す場合は入力個数を先に指定する必要がある。
    例
    ```console  
    # パイプ
    $ cat accession.txt | python fetch_gbk.py --email ~ --stdin

    # 非パイプ
    $ python fetch_gbk.py --email ~ --stdin
    1
    NC_012920
    ```
- `accessions`
    取得したいデータの開くセッション番号（`--stdin`オプションを使わない場合に1個以上指定する必要がある。）

## 動作条件
- python3.x
- biopython
