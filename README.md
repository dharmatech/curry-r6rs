# curry-r6rs

Portable R6RS versions of libraries from [Curry Scheme](https://github.com/deconstructo/curry) that Curry implements in C.

The question this project explores: does Curry's math code need to be C? If a plain R6RS library runs as fast on a good compiler like Chez Scheme, it can be shared with every R6RS implementation instead of living inside one.

## Libraries

| Library | Ported from | Status |
|---|---|---|
| `(curry-r6rs multivector)` | `src/multivec.c` | done, benchmarked |

### `(curry-r6rs multivector)`

Multivectors in the Clifford algebra Cl(p,q,r): p generators square to +1, q to −1, r to 0. That covers complex numbers (Cl(0,1,0)), quaternions as rotors (Cl(3,0,0)), spacetime algebra (Cl(1,3,0)), projective geometric algebra (Cl(3,0,1)) and conformal geometric algebra (Cl(4,1,0)), up to 16 generators.

```scheme
(import (rnrs) (curry-r6rs multivector))

(define g3  (make-algebra 3 0 0))
(define e1  (mv-basis g3 1))
(define e2  (mv-basis g3 2))
(define e12 (mv* e1 e2))                ; e1e2

(mv-scalar (mv* e12 e12))               ; => -1.0
(multivector->list (mv-wedge e1 e2))    ; => (0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0)
```

Blades are numbered by bitmap: bit k set means generator e(k+1) is present, so the coefficients run scalar, e1, e2, e12, e3, e13, e23, e123, ...

| Procedure | |
|---|---|
| `(make-algebra p q r)` | the algebra Cl(p,q,r) |
| `algebra?` `algebra-p` `algebra-q` `algebra-r` | |
| `(algebra-generators alg)` | n = p+q+r |
| `(algebra-blade-count alg)` | 2^n |
| `(make-multivector alg)` | zero |
| `(mv-basis alg i ...)` | unit blade e_i e_j ...; indices 1-based, any order (`(mv-basis g3 2 1)` is −e12); no indices gives 1 |
| `(list->multivector alg coefficients)` | one real per blade, in bitmap order |
| `(multivector->list m)` | |
| `multivector?` `multivector-algebra` `multivector-copy` | |
| `(multivector-ref m blade)` `(multivector-set! m blade x)` | blade is a bitmap or a list of 1-based indices |
| `(mv+ a b ...)` `(mv- a b ...)` `(mv- a)` | |
| `(mv* a b ...)` | geometric product, left-associative |
| `(mv-scale m s)` | |
| `(mv-wedge a b)` | outer product |
| `(mv-lcontract a b)` | left contraction |
| `mv-reverse` `mv-involute` `mv-conjugate` | |
| `(mv-grade m k)` | grade-k part |
| `mv-scalar` `mv-norm2` `mv-norm` `mv-normalize` | |
| `(mv-dual m)` | m times the reverse of the pseudoscalar |

Differences from Curry's built-in version:

- An algebra is a value you create once, instead of passing `p q r` to every constructor.
- Up to 16 generators; Curry allows 8.
- `list->multivector` requires exactly one coefficient per blade; Curry ignores extra or missing ones.
- Coefficients must be real. Curry silently drops the imaginary part of a complex number.
- `(mv+ a)` and `(mv* a)` return a copy, so the result never aliases a mutable argument.
- No quaternion conversions, since R6RS has no quaternion type.
- No custom printer, since R6RS has no portable way to define one. Use `multivector->list`.

Otherwise the products do the same floating-point arithmetic as Curry's C code, in the same order. Only the integer work of finding each blade's sign differs: the C code moves generators one at a time, while this library uses a closed form with no data-dependent branches (see `blade-product`). `tests/differential.sh` runs every operation on both and compares the results exactly, including the sign of zero. They match bit for bit against a Curry built without FMA contraction, which is GCC's x86-64 default. Where the C compiler fuses multiply-adds (arm64, or `-march=native`), the last bit of some results differs. NaN sign and payload also may differ, since IEEE 754 leaves them unspecified.

## Results

[docs/performance/2026-11-12-multivector.md](docs/performance/2026-11-12-multivector.md) is a self-contained write-up of the first table below: method, workloads, and how the gap to C was narrowed.

Measured 2026-11-12 on an Intel Core i9-10885H running Linux, at commit `869b9c6`. Chez Scheme 10.4.1 ran `(curry-r6rs multivector)`; Curry 1.30.1 (Release build, GCC `-O2`, no JIT) ran its built-in C multivectors. Each figure is the best of 5 runs, per operation; repeated runs vary by about 5%. The raw output is in [`bench/results/`](bench/results/). All three runs print identical checksums, and `tests/differential.sh` confirms they compute the same bits.

| Workload | Curry (C) | Chez | Chez, optimize-level 3 |
|---|---|---|---|
| Geometric product, Cl(3,0,0), 8 coefficients | 0.70 µs | 1.52 µs (2.2× slower) | 0.84 µs (1.2× slower) |
| Geometric product, Cl(4,1,0), 32 coefficients | 10.0 µs | 30.0 µs (3.0× slower) | 17.2 µs (1.7× slower) |
| Geometric product, Cl(8,0,0), 256 coefficients | 1.20 ms | 2.40 ms (2.0× slower) | 1.58 ms (1.3× slower) |
| Rotate a vector by a rotor, Cl(3,0,0) | 1.27 µs | 0.92 µs (1.4× faster) | 0.52 µs (2.4× faster) |

What the numbers say:

- **On a single large product, C is still faster**: 1.3–1.7× at optimize-level 3, 2–3× in Chez's default safe mode. A large product is almost all inner loop. Chez allocates nothing inside it, so the remaining gap is code generation. One specific piece was in Chez itself: R6RS's `fxbit-count` is always an out-of-line procedure call, while Chez's own `fxpopcount` compiles inline. A small patch to Chez closes most of the remaining gap; see below.
- **The rotor workload is what application code looks like**: build a vector, apply a rotor, read the result, half a million times. Here the portable library on Chez is faster than Curry's C code. Each product is small, so the time goes to the code around it, which Chez compiles to native code while Curry runs it in its bytecode VM with every float allocated on the heap.
- **The first version was a line-by-line port of Curry's loop**, at commit `d113422`; [its results](bench/results/2026-11-12-X1E-d113422.txt) are kept. That loop branched on effectively random sign bits, which the CPU mispredicts about half the time. The current closed form has no such branches and cut the large products by 40–50%: Cl(8,0,0) went from 4.3 to 2.4 ms in safe mode and from 3.3 to 1.6 ms at optimize-level 3.
- **C could use the same closed form.** In a quick standalone test, it helped C only with a hardware popcount instruction (`-mpopcnt`), which Curry's default x86-64 build doesn't use; without one, Curry's current loop is the faster choice. A multiplication table precomputed once per algebra would likewise speed up both.

### With `fxbit-count` inlined in Chez

Chez Scheme 10.4.1 compiles `fxbit-count` (R6RS's standard bit count) as a call into a library procedure. Its own non-standard `fxpopcount`, on the other hand, compiles to inline code, a single `POPCNT` instruction on x86-64. The library's blade product counts bits in its innermost loop.

A small patch to Chez inlines `fxbit-count` the same way. It is proposed upstream as [cisco/ChezScheme#1099](https://github.com/cisco/ChezScheme/pull/1099), open as of 2026-11-12. The table below measured an earlier version of the patch, applied to the 10.4.1 release. The pull request has a revised, branch-free version against Chez's development branch; there, the pull request reports Cl(8,0,0) products going from 1.56 ms to 1.09 ms at optimize-level 3. To compare like with like, both Chez builds below come from the same source with the same compiler, and Curry ran in the same sitting. [Raw results](bench/results/2026-11-12-X1E-5cf8b76-chez-fxbit-count.txt).

| Workload, optimize-level 3 | Curry (C) | Chez 10.4.1 | Chez 10.4.1 + patch |
|---|---|---|---|
| Geometric product, Cl(3,0,0) | 0.66 µs | 0.84 µs | 0.63 µs (on par with C) |
| Geometric product, Cl(4,1,0) | 9.3 µs | 17.2 µs | 11.5 µs (1.2× slower than C) |
| Geometric product, Cl(8,0,0) | 1.15 ms | 1.57 ms | 1.02 ms (1.1× faster than C) |
| Rotate a vector by a rotor | 1.18 µs | 0.51 µs | 0.42 µs (2.8× faster than C) |

In Chez's default safe mode the patch helps too: 10–22% on the products and 9% on the rotor workload. With it, portable R6RS code on Chez is within about 20% of Curry's hand-written C on every raw product, and matches or beats it on two of the three.

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
