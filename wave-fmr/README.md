# wave-fmr
Support schemes:

- FTCS
- Lax
- Lax-Wendroff
- MacCormack
- Upwind

## あそびかた

```
$ docker run -it --rm -v $PWD:/work -u $UID:$GID ishiy1993/formura-bin ./run upwind 0.5
$ ./plot data/upwind-0.5.dat
$ evince data/upwind-0.5.pdf
```

# Memo
wave-hsではなかった端の振動はなぜ生じるのか?
