# curry-r6rs

Portable R6RS versions of libraries from [Curry Scheme](https://github.com/deconstructo/curry) that Curry implements in C.

The question this project explores: does Curry's math code need to be C? If a plain R6RS library runs as fast on a good compiler like Chez Scheme, it can be shared with every R6RS implementation instead of living inside one.

## Libraries

| Library | What it provides | Ported from | Status |
|---|---|---|---|
| [`(curry-r6rs multivector)`](docs/multivector.md) | Clifford algebra multivectors, Cl(p,q,r) | [`src/multivec.c`](https://github.com/deconstructo/curry/blob/main/src/multivec.c) | done, benchmarked |

## Performance

On Chez Scheme 10.4.1 at optimize-level 3, `(curry-r6rs multivector)` is 2.4× faster than Curry's C code at rotating vectors, the kind of work application code does. On a single large geometric product, it is 1.2–1.7× slower. Both compute bit-identical results.

- [Multivector library vs Curry's C](docs/performance/2026-11-12-multivector.md) (2026-11-12): method, results, and how the gap to C was narrowed.
- [With `fxbit-count` inlined in Chez](docs/performance/2026-11-12-chez-fxbit-count.md) (2026-11-12): a proposed Chez patch, [cisco/ChezScheme#1099](https://github.com/cisco/ChezScheme/pull/1099), that brings every product to within 25% of C.

## Running

Requires [Chez Scheme](https://cisco.github.io/ChezScheme/) 10. From the repository root:

```sh
scheme --program tests/multivector.sps
```

These need a built Curry, by default at `~/src/curry/build-release/curry` (set `CURRY` to use another):

```sh
tests/differential.sh   # compare results with Curry's C code
bench/run.sh            # benchmark both; results go to bench/results/
```

## License

GPL-3.0-or-later, the same as Curry, since this is a port of Curry's code.
