# `(curry-r6rs multivector)`

A port of Curry Scheme's [`src/multivec.c`](https://github.com/deconstructo/curry/blob/main/src/multivec.c).

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

## Blades

Blades are numbered by bitmap: bit k set means generator e(k+1) is present, so the coefficients run scalar, e1, e2, e12, e3, e13, e23, e123, ...

## Procedures

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

## Differences from Curry's built-in version

- An algebra is a value you create once, instead of passing `p q r` to every constructor.
- Up to 16 generators; Curry allows 8.
- `list->multivector` requires exactly one coefficient per blade; Curry ignores extra or missing ones.
- Coefficients must be real. Curry silently drops the imaginary part of a complex number.
- `(mv+ a)` and `(mv* a)` return a copy, so the result never aliases a mutable argument.
- No quaternion conversions, since R6RS has no quaternion type.
- No custom printer, since R6RS has no portable way to define one. Use `multivector->list`.

## Agreement with Curry's results

The products do the same floating-point arithmetic as Curry's C code, in the same order. Only the integer work of finding each blade's sign differs: the C code moves generators one at a time, while this library uses a closed form with no data-dependent branches (see `blade-product` in [`multivector.sls`](../curry-r6rs/multivector.sls)). [`tests/differential.sh`](../tests/differential.sh) runs every operation on both and compares the results exactly, including the sign of zero. They match bit for bit against a Curry built without FMA contraction, which is GCC's x86-64 default. Where the C compiler fuses multiply-adds (arm64, or `-march=native`), the last bit of some results differs. NaN sign and payload also may differ, since IEEE 754 leaves them unspecified.

## Performance

See the [performance write-ups](../README.md#performance).
