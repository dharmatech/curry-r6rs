#!r6rs
;;; (curry-r6rs multivector) -- multivectors in the Clifford algebra Cl(p,q,r)
;;;
;;; A port of Curry Scheme's src/multivec.c to portable R6RS.
;;;   Curry Scheme: Copyright (C) 2026 Scáth
;;;   This port:    Copyright (C) 2026 Eduardo Cavazos
;;;
;;; This program is free software: you can redistribute it and/or modify it
;;; under the terms of the GNU General Public License as published by the Free
;;; Software Foundation, either version 3 of the License, or (at your option)
;;; any later version.  See the LICENSE file for details.
;;;
;;; ---------------------------------------------------------------------------
;;;
;;; An algebra Cl(p,q,r) has n = p+q+r generators e1 ... en:
;;;
;;;   the first p square to +1, the next q square to -1, the last r square to 0.
;;;
;;; A multivector holds one flonum coefficient per basis blade, 2^n in all.
;;; Blades are numbered by bitmap: bit k set means generator e(k+1) is present,
;;; so blade 0 is the scalar, 1 is e1, 2 is e2, 3 is e1e2 (e12), and so on.
;;; Coefficients live in a bytevector as native IEEE doubles, so they are
;;; stored unboxed.
;;;
;;; The products use the same algorithm as Curry's C code, including skipping
;;; zero coefficients, so for the same inputs the results agree bit for bit.

(library (curry-r6rs multivector)
  (export
    ;; algebras
    make-algebra algebra? algebra-p algebra-q algebra-r
    algebra-generators algebra-blade-count
    ;; construction and access
    make-multivector multivector? multivector-algebra
    multivector-ref multivector-set! multivector-copy
    mv-basis list->multivector multivector->list
    ;; arithmetic
    mv+ mv- mv* mv-scale mv-wedge mv-lcontract
    ;; involutions
    mv-reverse mv-involute mv-conjugate
    ;; projections and norms
    mv-dual mv-grade mv-scalar mv-norm2 mv-norm mv-normalize)

  (import (rnrs))

  ;; 2^16 coefficients is 512 KiB per multivector, and the geometric product
  ;; is O(4^n).  The limit also keeps every blade bitmap a small fixnum.
  (define max-generators 16)

  ;; -------------------------------------------------------------------------
  ;; Algebras

  (define-record-type algebra
    (fields p q r generators blade-count)
    (protocol
      (lambda (new)
        (lambda (p q r)
          (define (count? x)
            (and (fixnum? x) (fx<=? 0 x max-generators)))
          (unless (and (count? p) (count? q) (count? r)
                       (fx<=? (fx+ p (fx+ q r)) max-generators))
            (assertion-violation 'make-algebra
              "signature must be non-negative with p+q+r <= 16" p q r))
          (let ((n (fx+ p (fx+ q r))))
            (new p q r n (fxarithmetic-shift-left 1 n))))))
    (sealed #t)
    (opaque #t)
    (nongenerative curry-r6rs:algebra))

  ;; Two algebras with the same signature are interchangeable.
  (define (same-algebra? a b)
    (or (eq? a b)
        (and (fx=? (algebra-p a) (algebra-p b))
             (fx=? (algebra-q a) (algebra-q b))
             (fx=? (algebra-r a) (algebra-r b)))))

  ;; -------------------------------------------------------------------------
  ;; Multivectors

  (define-record-type (multivector %make-multivector multivector?)
    (fields (immutable algebra multivector-algebra)
            (immutable coefficients coefficients))
    (sealed #t)
    (opaque #t)
    (nongenerative curry-r6rs:multivector))

  (define-syntax cref
    (syntax-rules ()
      ((_ bv i) (bytevector-ieee-double-native-ref bv (fx* 8 i)))))

  (define-syntax cset!
    (syntax-rules ()
      ((_ bv i x) (bytevector-ieee-double-native-set! bv (fx* 8 i) x))))

  ;; All-zero bytes are +0.0 in IEEE 754.
  (define (zero-like alg)
    (%make-multivector alg (make-bytevector (fx* 8 (algebra-blade-count alg)) 0)))

  (define (check-algebra who x)
    (unless (algebra? x)
      (assertion-violation who "not an algebra" x)))

  (define (check-mv who x)
    (unless (multivector? x)
      (assertion-violation who "not a multivector" x)))

  (define (check-same who a b)
    (check-mv who a)
    (check-mv who b)
    (unless (same-algebra? (multivector-algebra a) (multivector-algebra b))
      (assertion-violation who "multivectors belong to different algebras" a b)))

  (define (check-real who x)
    (unless (real? x)
      (assertion-violation who "not a real number" x)))

  (define (make-multivector alg)
    (check-algebra 'make-multivector alg)
    (zero-like alg))

  ;; A blade is named either by its bitmap or by a list of 1-based generator
  ;; indices.  A list names the blade as a set: '(2 1) and '(1 2) both mean
  ;; blade 3, and the order does not change any sign.
  (define (blade-index who alg spec)
    (cond
      ((fixnum? spec)
       (unless (and (fx>=? spec 0) (fx<? spec (algebra-blade-count alg)))
         (assertion-violation who "blade index out of range" spec))
       spec)
      ((list? spec)
       (let ((n (algebra-generators alg)))
         (let loop ((s spec) (bits 0))
           (if (null? s)
               bits
               (let ((k (car s)))
                 (unless (and (fixnum? k) (fx<=? 1 k n))
                   (assertion-violation who "generator index out of range" k spec))
                 (when (fxbit-set? bits (fx- k 1))
                   (assertion-violation who "repeated generator index" k spec))
                 (loop (cdr s)
                       (fxior bits (fxarithmetic-shift-left 1 (fx- k 1)))))))))
      (else
       (assertion-violation who "blade must be a bitmap or a list of indices" spec))))

  (define (multivector-ref m blade)
    (check-mv 'multivector-ref m)
    (cref (coefficients m)
          (blade-index 'multivector-ref (multivector-algebra m) blade)))

  (define (multivector-set! m blade x)
    (check-mv 'multivector-set! m)
    (check-real 'multivector-set! x)
    (cset! (coefficients m)
           (blade-index 'multivector-set! (multivector-algebra m) blade)
           (real->flonum x)))

  (define (multivector-copy m)
    (check-mv 'multivector-copy m)
    (%make-multivector (multivector-algebra m) (bytevector-copy (coefficients m))))

  ;; Number of out-of-order pairs in a list of fixnums.  Its parity is the
  ;; sign picked up by sorting the generators into canonical order.
  (define (inversions lst)
    (let outer ((l lst) (count 0))
      (if (null? l)
          count
          (outer (cdr l)
                 (let inner ((r (cdr l)) (c count))
                   (cond ((null? r) c)
                         ((fx>? (car l) (car r)) (inner (cdr r) (fx+ c 1)))
                         (else (inner (cdr r) c))))))))

  ;; (mv-basis alg i ...) is the unit blade e_i e_j ...  With no indices it is
  ;; the scalar 1.  Indices are 1-based and may come in any order; (mv-basis
  ;; alg 2 1) is e2e1 = -e12.
  (define (mv-basis alg . indices)
    (check-algebra 'mv-basis alg)
    (let ((bits (blade-index 'mv-basis alg indices))
          (m (zero-like alg)))
      (cset! (coefficients m) bits (if (fxeven? (inversions indices)) 1.0 -1.0))
      m))

  ;; Coefficients in blade-bitmap order: scalar, e1, e2, e12, e3, ...
  (define (list->multivector alg lst)
    (check-algebra 'list->multivector alg)
    (let ((size (algebra-blade-count alg)))
      (unless (and (list? lst) (= (length lst) size))
        (assertion-violation 'list->multivector
          "expected a list with one coefficient per blade" size lst))
      (let* ((m (zero-like alg))
             (bv (coefficients m)))
        (let loop ((l lst) (i 0))
          (unless (null? l)
            (check-real 'list->multivector (car l))
            (cset! bv i (real->flonum (car l)))
            (loop (cdr l) (fx+ i 1))))
        m)))

  (define (multivector->list m)
    (check-mv 'multivector->list m)
    (let ((bv (coefficients m)))
      (let loop ((i (fx- (algebra-blade-count (multivector-algebra m)) 1))
                 (acc '()))
        (if (fx<? i 0)
            acc
            (loop (fx- i 1) (cons (cref bv i) acc))))))

  ;; -------------------------------------------------------------------------
  ;; Element-wise arithmetic

  (define-syntax define-elementwise
    (syntax-rules ()
      ((_ name who op)
       (define (name a b)
         (check-same who a b)
         (let* ((alg (multivector-algebra a))
                (size (algebra-blade-count alg))
                (ac (coefficients a))
                (bc (coefficients b))
                (res (zero-like alg))
                (rc (coefficients res)))
           (do ((i 0 (fx+ i 1)))
               ((fx=? i size) res)
             (cset! rc i (op (cref ac i) (cref bc i)))))))))

  (define-elementwise add2 'mv+ fl+)
  (define-elementwise sub2 'mv- fl-)

  (define (negate a)
    (check-mv 'mv- a)
    (let* ((alg (multivector-algebra a))
           (size (algebra-blade-count alg))
           (ac (coefficients a))
           (res (zero-like alg))
           (rc (coefficients res)))
      (do ((i 0 (fx+ i 1)))
          ((fx=? i size) res)
        (cset! rc i (fl- (cref ac i))))))

  (define (mv-scale a s)
    (check-mv 'mv-scale a)
    (check-real 'mv-scale s)
    (let* ((s (real->flonum s))
           (alg (multivector-algebra a))
           (size (algebra-blade-count alg))
           (ac (coefficients a))
           (res (zero-like alg))
           (rc (coefficients res)))
      (do ((i 0 (fx+ i 1)))
          ((fx=? i size) res)
        (cset! rc i (fl* (cref ac i) s)))))

  ;; -------------------------------------------------------------------------
  ;; Products

  ;; The product of basis blades a and b (bitmaps) in an algebra whose first
  ;; p generators square to +1 and next q to -1; the rest square to 0.
  ;;
  ;; Returns 0 when the product vanishes, otherwise +(c+1) or -(c+1) where c
  ;; is the result blade's bitmap and the sign is the product's sign.  One
  ;; fixnum, so the inner loops allocate nothing.
  ;;
  ;; Moves each generator k of b, lowest first, leftward into a: every
  ;; generator of a above k is one swap, one sign flip.  If a already holds
  ;; k, the pair squares to +1, -1 or 0 by k's position in the signature.
  (define (blade-product a b p q)
    (let ((pq (fx+ p q)))
      (let loop ((acc a) (b b) (neg? #f))
        (if (fx=? b 0)
            (if neg? (fx- -1 acc) (fx+ acc 1))
            (let* ((k (fxfirst-bit-set b))
                   (b (fxand b (fx- b 1)))
                   (neg? (if (fxodd? (fxbit-count (fxarithmetic-shift-right acc (fx+ k 1))))
                             (not neg?)
                             neg?)))
              (cond
                ((not (fxbit-set? acc k))
                 (loop (fxior acc (fxarithmetic-shift-left 1 k)) b neg?))
                ((fx>=? k pq) 0)
                (else
                 (loop (fxxor acc (fxarithmetic-shift-left 1 k))
                       b
                       (if (fx>=? k p) (not neg?) neg?)))))))))

  ;; Every product sums coefficient pairs through blade-product; they differ
  ;; only in which blade pairs (i, j) they skip.
  (define-syntax define-product
    (syntax-rules ()
      ((_ name who skip?)
       (define (name a b)
         (check-same who a b)
         (let* ((alg (multivector-algebra a))
                (p (algebra-p alg))
                (q (algebra-q alg))
                (size (algebra-blade-count alg))
                (ac (coefficients a))
                (bc (coefficients b))
                (res (zero-like alg))
                (rc (coefficients res)))
           (do ((i 0 (fx+ i 1)))
               ((fx=? i size) res)
             (let ((x (cref ac i)))
               (unless (fl=? x 0.0)
                 (do ((j 0 (fx+ j 1)))
                     ((fx=? j size))
                   (let ((y (cref bc j)))
                     (unless (or (fl=? y 0.0) (skip? i j))
                       (let ((s (blade-product i j p q)))
                         (cond
                           ((fx>? s 0)
                            (let ((k (fx- s 1)))
                              (cset! rc k (fl+ (cref rc k) (fl* x y)))))
                           ((fx<? s 0)
                            (let ((k (fx- -1 s)))
                              (cset! rc k (fl- (cref rc k) (fl* x y))))))))))))))))))

  ;; Geometric product: every pair of blades.
  (define-product geometric2 'mv* (lambda (i j) #f))

  ;; Outer product: blades sharing a generator contribute nothing.
  (define-product mv-wedge 'mv-wedge
    (lambda (i j) (not (fx=? (fxand i j) 0))))

  ;; Left contraction: only pairs where blade i is contained in blade j.
  (define-product mv-lcontract 'mv-lcontract
    (lambda (i j) (not (fx=? (fxand i j) i))))

  ;; One argument returns a fresh copy, so the result never aliases an input.
  (define mv+
    (case-lambda
      ((a) (multivector-copy a))
      ((a b) (add2 a b))
      ((a b . rest) (fold-left add2 (add2 a b) rest))))

  ;; One argument negates; more subtract left to right.
  (define mv-
    (case-lambda
      ((a) (negate a))
      ((a b) (sub2 a b))
      ((a b . rest) (fold-left sub2 (sub2 a b) rest))))

  ;; Left-associative: (mv* a b c) is (a b) c.
  (define mv*
    (case-lambda
      ((a) (multivector-copy a))
      ((a b) (geometric2 a b))
      ((a b . rest) (fold-left geometric2 (geometric2 a b) rest))))

  ;; -------------------------------------------------------------------------
  ;; Involutions: each flips the sign of some grades.

  (define-syntax define-grade-sign
    (syntax-rules ()
      ((_ name who flip?)
       (define (name a)
         (check-mv who a)
         (let* ((alg (multivector-algebra a))
                (size (algebra-blade-count alg))
                (ac (coefficients a))
                (res (zero-like alg))
                (rc (coefficients res)))
           (do ((i 0 (fx+ i 1)))
               ((fx=? i size) res)
             (let ((x (cref ac i)))
               (cset! rc i (if (flip? (fxbit-count i)) (fl- x) x)))))))))

  ;; Reverse: grade k flips when k(k-1)/2 is odd (grades 2, 3, 6, 7, ...).
  (define-grade-sign mv-reverse 'mv-reverse
    (lambda (k) (fxodd? (fxarithmetic-shift-right (fx* k (fx- k 1)) 1))))

  ;; Grade involution: odd grades flip.
  (define-grade-sign mv-involute 'mv-involute
    (lambda (k) (fxodd? k)))

  ;; Clifford conjugate, the involute of the reverse: grade k flips when
  ;; k(k+1)/2 is odd (grades 1, 2, 5, 6, ...).
  (define-grade-sign mv-conjugate 'mv-conjugate
    (lambda (k) (fxodd? (fxarithmetic-shift-right (fx* k (fx+ k 1)) 1))))

  ;; -------------------------------------------------------------------------
  ;; Projections and norms

  ;; Right complement: a times the reverse of the pseudoscalar e1...en.  In a
  ;; degenerate algebra (r > 0) the pseudoscalar squares to 0, so this is the
  ;; metric dual and not the Poincare/Hodge dual that PGA usually wants.
  (define (mv-dual a)
    (check-mv 'mv-dual a)
    (let* ((alg (multivector-algebra a))
           (pseudo (zero-like alg)))
      (cset! (coefficients pseudo) (fx- (algebra-blade-count alg) 1) 1.0)
      (geometric2 a (mv-reverse pseudo))))

  ;; The grade-k part of a.  A k outside 0..n gives zero.
  (define (mv-grade a k)
    (check-mv 'mv-grade a)
    (unless (fixnum? k)
      (assertion-violation 'mv-grade "grade must be a fixnum" k))
    (let* ((alg (multivector-algebra a))
           (size (algebra-blade-count alg))
           (ac (coefficients a))
           (res (zero-like alg))
           (rc (coefficients res)))
      (do ((i 0 (fx+ i 1)))
          ((fx=? i size) res)
        (when (fx=? (fxbit-count i) k)
          (cset! rc i (cref ac i))))))

  (define (mv-scalar a)
    (check-mv 'mv-scalar a)
    (cref (coefficients a) 0))

  ;; The scalar part of reverse(a) * a.  Negative for some multivectors in
  ;; mixed-signature algebras.
  (define (mv-norm2 a)
    (check-mv 'mv-norm2 a)
    (cref (coefficients (geometric2 (mv-reverse a) a)) 0))

  (define (mv-norm a)
    (flsqrt (flabs (mv-norm2 a))))

  (define (mv-normalize a)
    (let ((n2 (mv-norm2 a)))
      (when (fl=? n2 0.0)
        (assertion-violation 'mv-normalize "multivector has zero norm" a))
      (mv-scale a (fl/ 1.0 (flsqrt (flabs n2))))))
)
