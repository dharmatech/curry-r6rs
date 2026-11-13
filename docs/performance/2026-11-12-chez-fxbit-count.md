# Multivector performance with `fxbit-count` inlined in Chez Scheme

*2026-11-12. Measured at commit `5cf8b76`. [Raw results](../../bench/results/2026-11-12-X1E-5cf8b76-chez-fxbit-count.txt).*

This follows [the multivector write-up](2026-11-12-multivector.md) of the same day, with the same library, workloads and machine. The one change is to Chez Scheme itself. Like that write-up, this is a snapshot as of its date.

## Background

The library's blade product counts bits in its innermost loop, using `fxbit-count`, R6RS's standard bit count. Chez Scheme 10.4.1 compiles every `fxbit-count` call as a call into a library procedure. Chez also has its own non-standard `fxpopcount`, which compiles to inline code: a single `POPCNT` instruction on x86-64. Portable R6RS code can't use `fxpopcount`.

## The patch

A small patch to Chez inlines `fxbit-count` the way `fxpopcount` is. It is proposed upstream as [cisco/ChezScheme#1099](https://github.com/cisco/ChezScheme/pull/1099), open as of 2026-11-12.

The numbers below measured an earlier version of the patch, applied to the 10.4.1 release. That version branched on the sign of the argument. The pull request has a revised, branch-free version against Chez's development branch; there, the pull request reports Cl(8,0,0) products going from 1.56 ms to 1.09 ms at optimize-level 3.

To compare like with like, both Chez builds were compiled from the same source with the same compiler, and Curry ran in the same sitting.

## Results

Time per operation at optimize-level 3:

| Workload | Curry (C) | Chez 10.4.1 | Chez 10.4.1 + patch |
|---|---|---|---|
| Geometric product, Cl(3,0,0) | 0.66 µs | 0.84 µs | 0.63 µs (on par with C) |
| Geometric product, Cl(4,1,0) | 9.3 µs | 17.2 µs | 11.5 µs (1.2× slower than C) |
| Geometric product, Cl(8,0,0) | 1.15 ms | 1.57 ms | 1.02 ms (1.1× faster than C) |
| Rotate a vector, Cl(3,0,0) | 1.18 µs | 0.51 µs | 0.42 µs (2.8× faster than C) |

At Chez's default optimize level:

| Workload | Curry (C) | Chez 10.4.1 | Chez 10.4.1 + patch |
|---|---|---|---|
| Geometric product, Cl(3,0,0) | 0.66 µs | 1.52 µs | 1.37 µs (−10%) |
| Geometric product, Cl(4,1,0) | 9.3 µs | 29.7 µs | 24.9 µs (−16%) |
| Geometric product, Cl(8,0,0) | 1.15 ms | 2.38 ms | 1.85 ms (−23%) |
| Rotate a vector, Cl(3,0,0) | 1.18 µs | 0.92 µs | 0.84 µs (−9%) |

At optimize-level 3, the patch cuts the products' time by 24–35%. With it, the portable R6RS library is within 25% of Curry's hand-written C on every product, and matches or beats C on two of the three. At the default level it helps less: 10–23% on the products.
