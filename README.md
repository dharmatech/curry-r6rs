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

Otherwise the products use the same algorithm as Curry's C code. `tests/differential.sh` runs every operation on both and compares the results exactly, including the sign of zero. They match bit for bit against a Curry built without FMA contraction, which is GCC's x86-64 default. Where the C compiler fuses multiply-adds (arm64, or `-march=native`), the last bit of some results differs. NaN sign and payload also may differ, since IEEE 754 leaves them unspecified.

## Results

Measured 2026-11-12 on an Intel Core i9-10885H running Linux. Chez Scheme 10.4.1 ran `(curry-r6rs multivector)`; Curry 1.30.1 (Release build, GCC `-O2`, no JIT) ran its built-in C multivectors. Each figure is the best of 5 runs, per operation. The raw output is in [`bench/results/`](bench/results/). All three runs print identical checksums, and `tests/differential.sh` confirms they compute the same bits.

| Workload | Curry (C) | Chez | Chez, optimize-level 3 |
|---|---|---|---|
| Geometric product, Cl(3,0,0), 8 coefficients | 0.68 µs | 1.47 µs (2.2× slower) | 0.99 µs (1.5× slower) |
| Geometric product, Cl(4,1,0), 32 coefficients | 9.4 µs | 49 µs (5.2× slower) | 35 µs (3.7× slower) |
| Geometric product, Cl(8,0,0), 256 coefficients | 1.17 ms | 4.3 ms (3.7× slower) | 3.3 ms (2.8× slower) |
| Rotate a vector by a rotor, Cl(3,0,0) | 1.20 µs | 0.89 µs (1.35× faster) | 0.58 µs (2.1× faster) |

What the numbers say:

- **A large product is almost all inner loop, and there C is 2–5× faster** than the same algorithm on Chez. Chez isn't allocating inside that loop (a product allocates only its result), and inlining the blade product by hand changed nothing. The gap is in how Chez compiles the bit-twiddling and data-dependent sign branches compared with GCC.
- **The rotor workload is what application code looks like**: build a vector, apply a rotor, read the result, half a million times. Here the portable library on Chez is faster than Curry's C code. Each product is small, so the time goes to the code around it, which Chez compiles to native code while Curry runs it in its bytecode VM with every float allocated on the heap.
- **This is a direct port of the same algorithm.** A multiplication table precomputed once per algebra would remove most of the inner-loop cost on Chez, but the C code could use one too.

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
