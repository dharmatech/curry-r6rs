# Multivector performance: portable R6RS on Chez Scheme vs Curry's C

*2026-11-12. Measured at commit `869b9c6`. [Raw results](../../bench/results/2026-11-12-X1E-869b9c6.txt).*

This is a snapshot of the code and tools as of the date above. Later measurements will get a write-up of their own instead of replacing this one.

## Summary

[Curry Scheme](https://github.com/deconstructo/curry) implements Clifford algebra multivectors in C ([`src/multivec.c`](https://github.com/deconstructo/curry/blob/main/src/multivec.c)). `(curry-r6rs multivector)` ports that code to portable R6RS. It imports nothing but `(rnrs)`, and its products do the same floating-point arithmetic in the same order as the C code, so the two give bit-identical results.

On Chez Scheme 10.4.1, the R6RS library:

- **is faster than Curry's C on application-style code.** Rotating a vector by a rotor takes 0.92 µs, against 1.27 µs for Curry. At Chez's optimize-level 3 it takes 0.52 µs, 2.4× faster than Curry.
- **is slower on a single large geometric product**, where all the time goes to one tight loop: 2–3× slower than C at Chez's default optimize level, and 1.2–1.7× slower at optimize-level 3.

Both systems are their released versions, unchanged. Everything here comes from work on the Scheme library itself.

## What was compared

| | Curry | curry-r6rs |
|---|---|---|
| Implementation | Curry 1.30.1, Release build: GCC 15.2 at `-O2`, default x86-64 target, without the optional LLVM JIT | Chez Scheme 10.4.1 |
| Multivector code | built in, C | `(curry-r6rs multivector)`, portable R6RS |
| Coefficient storage | C array of `double` | bytevector of IEEE doubles (unboxed) |

**The same algorithm.** Both compute a product by visiting every pair of nonzero coefficients and adding ±a·b into the resulting blade. Neither uses a precomputed multiplication table, though either could. The one difference is how each pair's sign is found; see [How the gap was narrowed](#how-the-gap-was-narrowed).

**The same results.** [`tests/differential.sh`](../../tests/differential.sh) runs every operation on both implementations and compares the results exactly, including the sign of zero. It covers 12 algebras with three kinds of input each: dense, with zeros and negative zeros, and with infinity and NaN. That's 712 results and 44,666 numbers, every one bit-identical, with NaN counted equal to NaN. The benchmarks' checksums also match on both sides.

**Machine.** Intel Core i9-10885H (8 cores, 2.4 GHz base), Linux 7.0.

**Method.** Each workload runs 5 times, with a garbage collection before each run. The tables show the best run divided by the number of operations. Medians are within 3% of the best runs, except Curry's rotor run at 7%; the raw results include both. Curry, which ran unchanged in both sittings described here, differed between them by up to 6%.

**Optimize levels.** Chez's default optimize level keeps all of Chez's own run-time type checks, and it's what a library user gets without asking. Optimize-level 3 (`scheme --optimize-level 3`) drops those checks, trusting the program to pass the right types.

## Workloads

| Workload | What it does | Operations |
|---|---|---|
| Geometric product, Cl(3,0,0) | product of two fixed multivectors, all 8 coefficients nonzero: 64 coefficient pairs | 1,000,000 |
| Geometric product, Cl(4,1,0) | the same with 32 coefficients: 1,024 pairs | 200,000 |
| Geometric product, Cl(8,0,0) | the same with 256 coefficients: 65,536 pairs | 1,000 |
| Rotate a vector, Cl(3,0,0) | build a vector from a list, rotate it by a fixed rotor (R v R̃, two products), add up three of its coefficients | 500,000 |

The code is in [`bench/multivector.sps`](../../bench/multivector.sps) (Chez) and [`bench/multivector-curry-c.scm`](../../bench/multivector-curry-c.scm) (Curry). The two use the same inputs and the same order of operations.

## Results

Time per operation:

| Workload | Curry (C) | Chez | Chez, optimize-level 3 |
|---|---|---|---|
| Geometric product, Cl(3,0,0) | 0.70 µs | 1.52 µs (2.2× slower) | 0.84 µs (1.2× slower) |
| Geometric product, Cl(4,1,0) | 10.0 µs | 30.0 µs (3.0× slower) | 17.2 µs (1.7× slower) |
| Geometric product, Cl(8,0,0) | 1.20 ms | 2.40 ms (2.0× slower) | 1.58 ms (1.3× slower) |
| Rotate a vector, Cl(3,0,0) | 1.27 µs | 0.92 µs (1.4× faster) | 0.52 µs (2.4× faster) |

**On a single large product, C is ahead.** A Cl(8,0,0) product is 65,536 multiply-adds in one loop nest and very little else. Chez allocates nothing inside that loop except the result multivector, so the remaining gap is in the machine code each compiler generates for it.

**The rotor workload looks like application code.** Its products are small, so most of the time goes to the code around them: building a list, making a multivector from it, reading coefficients back out, and adding floats. Chez compiles all of that to native code. Curry runs it in its bytecode VM, where every float is allocated on the heap. A Curry built with its LLVM JIT might do better here; that wasn't measured.

## How the gap was narrowed

The first version of the library ported Curry's sign computation line by line. To multiply blade `a` by blade `b`, that loop moves `b`'s generators into place one at a time. For each generator, it counts the generators it passes, flips the sign when the count is odd, and checks whether the generator meets itself. The results of that version, at commit `d113422`, are [kept](../../bench/results/2026-11-12-X1E-d113422.txt). At optimize-level 3 it was 3.7× slower than C on Cl(4,1,0) and 2.8× slower on Cl(8,0,0).

The time was not going to allocation or procedure calls. It went to branch mispredictions. Whether a sign flips depends on the data, and across all the blade pairs of a product it is effectively random, so the CPU guesses wrong about half the time.

The current version computes the same answer in closed form, with no branch on the sign:

- The resulting blade is `a xor b`, since generators the two have in common cancel out.
- The product is zero if a common generator squares to 0.
- Sorting the generators of `a` followed by those of `b` into order takes one swap for every generator of `a` that is greater than a generator of `b`. That count is the sum of `popcount((a >> d) & b)` over d ≥ 1. Each common generator that squares to −1 adds one more sign flip. An odd total means a negative sign.

```scheme
;; Returns -1 if the product is zero, otherwise (blade << 1) | parity.
(define (blade-product a b p q)
  (let ((pq (fx+ p q))
        (shared (fxand a b)))
    (if (fx=? (fxarithmetic-shift-right shared pq) 0)
        (let loop ((x (fxarithmetic-shift-right a 1))
                   ;; start from the shared generators that square to -1
                   (count (fxbit-count
                            (fxand (fxarithmetic-shift-right shared p)
                                   (fx- (fxarithmetic-shift-left 1 q) 1)))))
          (if (fx=? x 0)
              (fxior (fxarithmetic-shift-left (fxxor a b) 1) (fxand count 1))
              (loop (fxarithmetic-shift-right x 1)
                    (fx+ count (fxbit-count (fxand x b))))))
        -1)))
```

The caller turns the parity into a multiplier, 1 − 2·parity, which is +1.0 or −1.0, and multiplies by it instead of branching. The loop here depends only on how many generators `a` has, not on the signs. This is about as many lines as the loop it replaced.

Only integer work changed. The floating-point operations, (±1.0 · x) · y added into the result, are the ones the C code performs, in the same order. The differential test still finds every result bit-identical.

Time per operation, before and after:

| Workload | Chez, before → after | Optimize-level 3, before → after |
|---|---|---|
| Geometric product, Cl(3,0,0) | 1.47 → 1.52 µs | 0.99 → 0.84 µs (−15%) |
| Geometric product, Cl(4,1,0) | 49.4 → 30.0 µs (−39%) | 35.4 → 17.2 µs (−51%) |
| Geometric product, Cl(8,0,0) | 4.29 → 2.40 ms (−44%) | 3.26 → 1.58 ms (−52%) |
| Rotate a vector, Cl(3,0,0) | 0.89 → 0.92 µs | 0.58 → 0.52 µs (−10%) |

The large products took 40–50% less time. The small ones changed little, since per-call costs dominate them. At Chez's default level they moved within the run-to-run noise.

## Would the C code benefit too?

The same closed form is available to C. A quick standalone test of Curry's inner loop (not part of this repository) gave mixed results. With `-mpopcnt`, which lets GCC use the x86-64 `POPCNT` instruction, the closed form made Cl(8,0,0) about 2.5× faster. Without it, GCC computes popcounts in software, and the closed form came out slower than Curry's current loop. GCC's default x86-64 target, which Curry's Release build uses, is the second case.

So Curry's loop is the right choice for its default build, and the comparison above is against Curry as it builds by default. A multiplication table precomputed once per algebra would speed up both sides. It was left out here to keep the comparison like for like.

## Reproducing

You need Chez Scheme 10 and a Curry Release build (`cmake -B build-release -DCMAKE_BUILD_TYPE=Release`). From the repository root, at commit `869b9c6`:

```sh
scheme --program tests/multivector.sps   # the library's own tests
tests/differential.sh                    # compare results with Curry's C code
bench/run.sh                             # benchmark both; results go to bench/results/
```

Both scripts look for Curry at `~/src/curry/build-release/curry`. Set `CURRY` to use another binary.
