# curry-r6rs

Portable R6RS versions of libraries from [Curry Scheme](https://github.com/deconstructo/curry) that Curry implements in C.

The question this project explores: does Curry's math code need to be C? If a plain R6RS library runs as fast on a good compiler like Chez Scheme, it can be shared with every R6RS implementation instead of living inside one.

## Libraries

| Library | Ported from | Status |
|---|---|---|
| `(curry-r6rs multivector)` | `src/multivec.c` | done |

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
- `list->multivector` requires exactly one coefficient per blade; Curry ignores extra or missing ones.
- `(mv+ a)` and `(mv* a)` return a copy, so the result never aliases a mutable argument.
- No quaternion conversions, since R6RS has no quaternion type.
- No custom printer, since R6RS has no portable way to define one. Use `multivector->list`.

Otherwise the products use the same algorithm as Curry's C code, so the same inputs give bit-identical results.

## Running

Requires [Chez Scheme](https://cisco.github.io/ChezScheme/) 10. From the repository root:

```sh
scheme --program tests/multivector.sps
```

## License

GPL-3.0-or-later, the same as Curry, since this is a port of Curry's code.
