# calculate_weights.py

## 使い方

```console
$ python calculate_coordinates.py gbkfiles... --weight weight.json --mapping mapping.json
```

- `--weight` 
    座標計算に使う重みを`json`ファイルで指定する。
    指定しない場合、重みは付与されない。
- `--mapping` 
    割り振るベクトルを`json`ファイルで指定する。
- `-d`, `--destination`
    指定した時、その名前のディレクトリに座標データを書き出す。
    指定されていない時、カレントディレクトリに`calculated`ディレクトリを作成する。

## 動作条件
- python3.x
- biopython
