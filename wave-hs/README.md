# wave-hs
Support schemes:

- FTCS
- Lax
- Lax-Wendroff
- Upwind
- MUSCL

```
$ stack exec -- wave-hs ftcs 0.5
$ ./plot "data/ftcs-0.5.dat"
$ evince "data/ftcs-0.5.pdf"
```
