# calculate_weights.py

## 使い方

```console
$ python calculate_weights.py gbkfiles... --allow_chars ATGC
```

- `--stdout`
    このフラグが指定された時、計算した重みを標準出力に書き出す
- `-d`, `--destination`
    このフラグら指定された時、そのディレクトリを作成し重みファイルを書き出す。
    指定されていない時、カレントディレクトリに`weights.json`で書き出す。
- `--allow_chars`
    重みを算出するときに許可する文字の集合。
    指定しない場合は`ATGC`の文字でフィルタされる。

## 動作条件
- python3.x
- biopython
