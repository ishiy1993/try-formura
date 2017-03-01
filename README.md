# Formuraであそぶ
[Formura](https://github.com/nushio3/formura)

## あそびかた
`hoge.fmr`と`hoge.yaml`を作成する。

```
$ docker run -it --rm -v $PWD:/src -w /src -u $UID:$GID ishiy1993/formura bash
$ /work/bin/formura hoge.fmr
$ mpicxx hoge-main.cpp hoge.c hoge_internal_*.c
$ mpirun -n 1 ./a.out
```

もしセグフォル場合は、出力ファイルを書き出すディレクトリが存在するか確認すること

