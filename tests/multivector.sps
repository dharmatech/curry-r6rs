#!r6rs
;;; Tests for (curry-r6rs multivector).
;;;
;;; The first section ports the multivector checks from Curry's
;;; tests/numeric_ext_tests.scm (all except printing and the quaternion
;;; round trip, which have no R6RS counterpart).  Two of Curry's checks
;;; compared only the scalar part, which is 0 whatever the sign; those now
;;; compare whole multivectors.  The rest are new.
;;;
;;; Run from the repository root:  scheme --program tests/multivector.sps

(import (rnrs) (curry-r6rs multivector))

;;; ---------------------------------------------------------------------------
;;; A minimal harness

(define pass 0)
(define fail 0)

(define (report ok? label detail)
  (if ok?
      (set! pass (+ pass 1))
      (begin
        (set! fail (+ fail 1))
        (display "FAIL: ") (display label)
        (when detail (display " -- ") (write detail))
        (newline))))

(define (check label actual expected)
  (report (equal? actual expected) label (list 'got actual 'expected expected)))

(define (check-approx label actual expected eps)
  (report (< (abs (- actual expected)) eps) label
          (list 'got actual 'expected expected)))

(define (check-true label x) (report x label #f))

;; Passes when thunk raises an assertion violation.
(define (check-raises label thunk)
  (report (guard (c ((assertion-violation? c) #t) (else #f))
            (thunk)
            #f)
          label
          "no assertion violation raised"))

;; Largest absolute difference between coefficients.
(define (distance a b)
  (fold-left (lambda (m x y) (max m (abs (- x y))))
             0.0
             (multivector->list a)
             (multivector->list b)))

(define (check-mv-approx label a b eps)
  (report (< (distance a b) eps) label
          (list (multivector->list a) (multivector->list b))))

;; Numeric equality per coefficient, so -0.0 matches 0.0.
(define (check-mv= label a b)
  (let ((la (multivector->list a)) (lb (multivector->list b)))
    (report (and (= (length la) (length lb)) (for-all = la lb))
            label (list 'got la 'expected lb))))

;; Exact per coefficient, including the sign of zero.
(define (check-eqv label actual expected)
  (report (eqv? actual expected) label (list 'got actual 'expected expected)))

;;; ---------------------------------------------------------------------------
;;; Ported from Curry: Cl(3,0,0)

(define g3   (make-algebra 3 0 0))
(define e1   (mv-basis g3 1))
(define e2   (mv-basis g3 2))
(define e3   (mv-basis g3 3))
(define e12  (mv-basis g3 1 2))
(define e23  (mv-basis g3 2 3))
(define e13  (mv-basis g3 1 3))
(define e123 (mv-basis g3 1 2 3))
(define one  (mv-basis g3))

(check "multivector? yes" (multivector? e1) #t)
(check "multivector? no"  (multivector? 42) #f)
(check "signature" (let ((a (multivector-algebra e1)))
                     (list (algebra-p a) (algebra-q a) (algebra-r a)))
       '(3 0 0))

(check "e1*e1 scalar" (mv-scalar (mv* e1 e1)) 1.0)
(check "e2*e2 scalar" (mv-scalar (mv* e2 e2)) 1.0)
(check "e3*e3 scalar" (mv-scalar (mv* e3 e3)) 1.0)

(define e12-prod (mv* e1 e2))
(define e21-prod (mv* e2 e1))
(check "e1*e2 grade-2 blade" (mv-scalar (mv* e12-prod e12-prod)) -1.0)
(check-mv= "e1*e2 = e12"         e12-prod e12)
(check-mv= "e1*e2 = -e2*e1"      e12-prod (mv- e21-prod))

(check-approx "mv-norm e1"  (mv-norm e1)  1.0 1e-12)
(check-approx "mv-norm e12" (mv-norm e12) 1.0 1e-12)
(check "mv-norm2 e1" (mv-norm2 e1) 1.0)

(check-approx "mv+ norm" (mv-norm (mv+ e1 e2)) (sqrt 2.0) 1e-12)
(check-approx "mv- recovers e1" (mv-norm (mv- (mv- (mv+ e1 e2) e2) e1)) 0.0 1e-12)
(check "mv-scale norm2" (mv-norm2 (mv-scale e1 2.0)) 4.0)

(define mixed (mv+ one e1 e12 e123))
(check "mv-grade 0" (mv-scalar (mv-grade mixed 0)) 1.0)
(check-approx "mv-grade 1 norm" (mv-norm (mv-grade mixed 1)) 1.0 1e-12)
(check-approx "mv-grade 2 norm" (mv-norm (mv-grade mixed 2)) 1.0 1e-12)
(check-approx "mv-grade 3 norm" (mv-norm (mv-grade mixed 3)) 1.0 1e-12)

(check-mv= "mv-reverse e12 = -e12" (mv-reverse e12) (mv- e12))
(check-approx "mv-involute e12"  (mv-norm (mv- (mv-involute e12) e12)) 0.0 1e-12)
(check-approx "mv-involute e1"   (mv-norm (mv+ (mv-involute e1) e1)) 0.0 1e-12)
(check-approx "mv-conjugate e12" (mv-norm (mv+ (mv-conjugate e12) e12)) 0.0 1e-12)

(check "wedge e1^e1 = 0"   (mv-scalar (mv-wedge e1 e1)) 0.0)
(check "wedge e1^e2 = e12" (mv-norm (mv- (mv-wedge e1 e2) e12)) 0.0)

(define from-list (list->multivector (make-algebra 2 0 0) '(1.0 2.0 3.0 4.0)))
(check "list->multivector multivector?" (multivector? from-list) #t)
(check "multivector-ref blade 0" (multivector-ref from-list 0) 1.0)
(check "multivector-ref blade 1" (multivector-ref from-list 1) 2.0)

(check-approx "mv-normalize" (mv-norm (mv-normalize (mv-scale e1 3.0))) 1.0 1e-12)

;;; ---------------------------------------------------------------------------
;;; Construction and access

(check "algebra-generators" (algebra-generators g3) 3)
(check "algebra-blade-count" (algebra-blade-count g3) 8)
(check "make-multivector is zero"
       (multivector->list (make-multivector g3)) '(0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0))
(check "mv-basis with no indices is 1" (multivector->list one)
       '(1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0))
(check-mv= "mv-basis 2 1 = -e12" (mv-basis g3 2 1) (mv- e12))
(check-mv= "mv-basis 2 1 3 = -e123" (mv-basis g3 2 1 3) (mv- e123))
(check-mv= "mv-basis 3 1 2 = e123" (mv-basis g3 3 1 2) e123)
(check "ref by index list" (multivector-ref e12 '(1 2)) 1.0)
(check "index list is a set" (multivector-ref e12 '(2 1)) 1.0)
(check "ref by empty list is the scalar" (multivector-ref one '()) 1.0)

(let ((m (make-multivector g3)))
  (multivector-set! m '(1 3) 5)
  (check "set! stores a flonum" (multivector-ref m 5) 5.0)
  (let ((c (multivector-copy m)))
    (multivector-set! c 5 7.0)
    (check "copy is independent" (multivector-ref m 5) 5.0)))

(check-true "unary mv+ copies" (not (eq? (mv+ e1) e1)))
(check-true "unary mv* copies" (not (eq? (mv* e1) e1)))
(check-mv= "unary mv- negates" (mv- e1) (mv-scale e1 -1))
(check-mv= "n-ary mv* is left-associative" (mv* e1 e2 e3) e123)
(check-mv= "n-ary mv- subtracts left to right" (mv- e1 e1 e1) (mv- e1))
(check-mv= "same signature, different algebra objects"
           (mv* e1 (mv-basis (make-algebra 3 0 0) 2)) e12)

;;; ---------------------------------------------------------------------------
;;; Products of basis blades in several signatures

(check-mv= "e1 e23 = e123" (mv* e1 e23) e123)
(check-mv= "e23 e1 = e123" (mv* e23 e1) e123)
(check-mv= "e12 e23 = e13" (mv* e12 e23) e13)
(check "e123^2 = -1 in Cl(3,0,0)" (mv-scalar (mv* e123 e123)) -1.0)

(let* ((c (make-algebra 0 1 0))
       (i (mv-basis c 1)))
  (check "e1^2 = -1 in Cl(0,1,0)" (multivector->list (mv* i i)) '(-1.0 0.0)))

(let* ((a (make-algebra 1 1 0)))
  (check "e1^2 = +1 in Cl(1,1,0)" (mv-scalar (mv* (mv-basis a 1) (mv-basis a 1))) 1.0)
  (check "e2^2 = -1 in Cl(1,1,0)" (mv-scalar (mv* (mv-basis a 2) (mv-basis a 2))) -1.0))

(let* ((pga (make-algebra 3 0 1))
       (e1 (mv-basis pga 1))
       (e4 (mv-basis pga 4)))
  (check "e4^2 = 0 in Cl(3,0,1)"
         (multivector->list (mv* e4 e4)) (multivector->list (make-multivector pga)))
  (check-mv= "e1 e4 = e14" (mv* e1 e4) (mv-basis pga 1 4))
  (check-mv= "e4 e1 = -e14" (mv* e4 e1) (mv- (mv-basis pga 1 4))))

(let* ((sta (make-algebra 1 3 0))
       (i (mv-basis sta 1 2 3 4)))
  (check "pseudoscalar^2 = -1 in Cl(1,3,0)" (mv-scalar (mv* i i)) -1.0))

;; Left contraction is zero unless the left blade is inside the right one.
(check-mv= "e1 contract e12 = e2"  (mv-lcontract e1 e12) e2)
(check-mv= "e2 contract e12 = -e1" (mv-lcontract e2 e12) (mv- e1))
(check-mv= "e12 contract e1 = 0"   (mv-lcontract e12 e1) (make-multivector g3))
(check-mv= "e1 contract e2 = 0"    (mv-lcontract e1 e2) (make-multivector g3))
(check-mv= "e1 contract e1 = 1"    (mv-lcontract e1 e1) one)
(check-mv= "1 contract e12 = e12"  (mv-lcontract one e12) e12)

;; The dual multiplies by the reverse of the pseudoscalar, -e123 here.
(check-mv= "dual 1 = -e123"   (mv-dual one) (mv- e123))
(check-mv= "dual e1 = -e23"   (mv-dual e1) (mv- e23))
(check-mv= "dual e12 = e3"    (mv-dual e12) e3)
(check-mv= "dual e123 = 1"    (mv-dual e123) one)
;; In Cl(3,0,1) the null generator e4 cancels against the pseudoscalar's e4.
(let ((pga (make-algebra 3 0 1)))
  (check-mv= "dual e1 = e234 in Cl(3,0,1)" (mv-dual (mv-basis pga 1)) (mv-basis pga 2 3 4))
  (check-mv= "dual e4 = 0 in Cl(3,0,1)" (mv-dual (mv-basis pga 4)) (make-multivector pga)))

;;; ---------------------------------------------------------------------------
;;; Zero-skipping and signed zeros, which must match Curry's C code exactly

;; (Not (mv-scale e1 +inf.0): that turns every zero into 0 * inf = NaN.)
(check-true "scaling by inf makes zeros NaN" (nan? (mv-scalar (mv-scale e1 +inf.0))))

(let* ((a (list->multivector g3 '(-0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0)))
       (b (list->multivector g3 '(0.0 +inf.0 0.0 0.0 0.0 0.0 0.0 0.0)))
       (p (mv* a b)))
  ;; -0.0 * inf would be NaN, but zero coefficients are skipped.
  (check-eqv "a -0.0 coefficient is skipped" (multivector-ref p 1) 0.0)
  (check-eqv "inf survives" (mv-scalar p) +inf.0))

(let ((p (mv* (list->multivector g3 '(+nan.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0)) e1)))
  (check-true "a NaN coefficient is not skipped" (nan? (multivector-ref p 1))))

(check-eqv "negating zero gives -0.0" (mv-scalar (mv- (make-multivector g3))) -0.0)
(check-eqv "reverse flips the sign of a zero bivector coefficient"
           (multivector-ref (mv-reverse (make-multivector g3)) 3) -0.0)
(check-eqv "reverse keeps the sign of a zero vector coefficient"
           (multivector-ref (mv-reverse (make-multivector g3)) 1) 0.0)

;;; ---------------------------------------------------------------------------
;;; The largest algebra

(let* ((g16 (make-algebra 16 0 0))
       (e16 (mv-basis g16 16)))
  (check "16 generators" (algebra-blade-count g16) 65536)
  (check "e16^2 = 1" (mv-scalar (mv* e16 e16)) 1.0)
  (check "e1 e16 = e1,16" (multivector-ref (mv* (mv-basis g16 1) e16) '(1 16)) 1.0)
  (check "e16 e1 = -e1,16" (multivector-ref (mv* e16 (mv-basis g16 1)) '(1 16)) -1.0))
(check "8+8 generators" (algebra-generators (make-algebra 8 8 0)) 16)

;;; ---------------------------------------------------------------------------
;;; Algebraic identities on pseudo-random dense multivectors

(define seed 12345)
(define (random-unit)                   ; a flonum in [-1, 1)
  (set! seed (mod (+ (* 1103515245 seed) 12345) 2147483648))
  (- (/ (inexact seed) 1073741824.0) 1.0))

(define (random-mv alg)
  (let loop ((i 0) (acc '()))
    (if (= i (algebra-blade-count alg))
        (list->multivector alg acc)
        (loop (+ i 1) (cons (random-unit) acc)))))

(define (random-vector alg)
  (mv-grade (random-mv alg) 1))

(define (iota* n)                       ; (0 1 ... n-1)
  (let loop ((i (- n 1)) (acc '()))
    (if (< i 0) acc (loop (- i 1) (cons i acc)))))

(for-each
  (lambda (sig)
    (let* ((alg (apply make-algebra sig))
           (name (lambda (s) (string-append s " in Cl" (call-with-string-output-port
                                                         (lambda (p) (write sig p))))))
           (a (random-mv alg)) (b (random-mv alg)) (c (random-mv alg))
           (u (random-vector alg)) (v (random-vector alg)))
      (check-mv-approx (name "associativity") (mv* (mv* a b) c) (mv* a (mv* b c)) 1e-9)
      (check-mv-approx (name "distributivity") (mv* a (mv+ b c)) (mv+ (mv* a b) (mv* a c)) 1e-9)
      (check-mv-approx (name "reverse is an anti-automorphism")
                       (mv-reverse (mv* a b)) (mv* (mv-reverse b) (mv-reverse a)) 1e-9)
      (check-mv-approx (name "involute is an automorphism")
                       (mv-involute (mv* a b)) (mv* (mv-involute a) (mv-involute b)) 1e-9)
      (check-mv= (name "conjugate = involute of reverse")
                 (mv-conjugate a) (mv-involute (mv-reverse a)))
      (check-mv-approx (name "wedge associativity")
                       (mv-wedge (mv-wedge a b) c) (mv-wedge a (mv-wedge b c)) 1e-9)
      (check-mv-approx (name "u v = u.v + u^v")
                       (mv* u v) (mv+ (mv-lcontract u v) (mv-wedge u v)) 1e-12)
      (check-mv-approx (name "vectors anticommute under wedge")
                       (mv-wedge u v) (mv- (mv-wedge v u)) 1e-12)
      (check-mv= (name "grades sum to the whole")
                 (fold-left mv+ (mv-grade a 0)
                            (map (lambda (k) (mv-grade a k))
                                 (cdr (iota* (+ 1 (algebra-generators alg))))))
                 a)))
  '((2 0 0) (3 0 0) (0 3 0) (1 3 0) (3 0 1) (4 1 0) (2 2 2)))

;;; ---------------------------------------------------------------------------
;;; Rotation with a rotor, after Curry's docs/reference/multivec.md

;; R = cos(h) - sin(h) B, where B = ax e23 + ay e31 + az e12 is the plane
;; dual to the axis.  e31 = -e13, so the e13 coefficient is +sin(h) ay.
;; (Curry's docs give it as -sin(h) ay, which turns the y axis the wrong way.)
(define (make-rotor ax ay az theta)
  (let ((h (/ theta 2.0))
        (m (make-multivector g3)))
    (multivector-set! m 0 (cos h))
    (multivector-set! m #b110 (* (- (sin h)) ax))
    (multivector-set! m #b101 (* (sin h) ay))
    (multivector-set! m #b011 (* (- (sin h)) az))
    m))

(define quarter-turn (* 2 (atan 1)))

;; Rotate vector v by a quarter turn about axis; return its e1, e2, e3 parts.
(define (quarter-turn-about axis v)
  (let* ((r (apply make-rotor (append axis (list quarter-turn))))
         (w (mv* r v (mv-reverse r))))
    (list (multivector-ref w 1) (multivector-ref w 2) (multivector-ref w 4))))

(define (check-rotation label actual expected)
  (report (for-all (lambda (x y) (< (abs (- x y)) 1e-12)) actual expected)
          label (list 'got actual 'expected expected)))

(check-rotation "about x, e2 turns to e3" (quarter-turn-about '(1 0 0) e2) '(0 0 1))
(check-rotation "about y, e3 turns to e1" (quarter-turn-about '(0 1 0) e3) '(1 0 0))
(check-rotation "about y, e1 turns to -e3" (quarter-turn-about '(0 1 0) e1) '(0 0 -1))
(check-rotation "about z, e1 turns to e2" (quarter-turn-about '(0 0 1) e1) '(0 1 0))
(check-approx "rotor is unit" (mv-norm2 (make-rotor 0.6 0 0.8 1.0)) 1.0 1e-12)

;;; ---------------------------------------------------------------------------
;;; Errors

(check-raises "make-algebra rejects negative" (lambda () (make-algebra -1 0 0)))
(check-raises "make-algebra rejects n > 16" (lambda () (make-algebra 9 8 0)))
(check-raises "make-algebra rejects non-fixnum" (lambda () (make-algebra 3.0 0 0)))
(check-raises "mixed algebras" (lambda () (mv* e1 (mv-basis (make-algebra 2 0 0) 1))))
(check-raises "blade bitmap out of range" (lambda () (multivector-ref e1 8)))
(check-raises "negative blade bitmap" (lambda () (multivector-ref e1 -1)))
(check-raises "generator index out of range" (lambda () (mv-basis g3 4)))
(check-raises "generator index zero" (lambda () (mv-basis g3 0)))
(check-raises "repeated generator index" (lambda () (mv-basis g3 1 1)))
(check-raises "list->multivector too short" (lambda () (list->multivector g3 '(1.0))))
(check-raises "list->multivector non-real" (lambda ()
                                             (list->multivector (make-algebra 1 0 0) '(1.0 x))))
(check-raises "set! non-real" (lambda () (multivector-set! (make-multivector g3) 0 'x)))
(check-raises "normalize zero" (lambda () (mv-normalize (make-multivector g3))))
(check-raises "not a multivector" (lambda () (mv* e1 42)))
(check-raises "mv-grade non-fixnum" (lambda () (mv-grade e1 1.0)))
(check-raises "make-multivector without algebra" (lambda () (make-multivector '(3 0 0))))
(check-raises "set! blade out of range" (lambda () (multivector-set! (make-multivector g3) 8 1.0)))
(check-raises "mv-scale non-real" (lambda () (mv-scale e1 'x)))
(check-raises "mv-scale complex" (lambda () (mv-scale e1 1+2i)))
(check-raises "unary mv+ non-multivector" (lambda () (mv+ 42)))
(let ((other (mv-basis (make-algebra 2 0 0) 1)))
  (for-each
    (lambda (name op)
      (check-raises (string-append name " with mixed algebras") (lambda () (op e1 other))))
    '("mv+" "mv-" "mv-wedge" "mv-lcontract")
    (list mv+ mv- mv-wedge mv-lcontract)))

;;; ---------------------------------------------------------------------------

(display pass) (display " passed, ") (display fail) (display " failed") (newline)
(exit (if (= fail 0) 0 1))
