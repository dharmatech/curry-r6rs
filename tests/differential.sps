#!r6rs
;;; Differential test, Chez side: print every operation's result as data.
;;;
;;; tests/differential-curry-c.scm prints the same thing using Curry's
;;; built-in C multivectors; tests/differential.sh compares the two exactly.
;;; Keep the two files in step: same inputs, same operations, same order.
;;;
;;; Each output line is (operation signature input-kind value ...).

(import (rnrs) (curry-r6rs multivector))

(define seed 20261112)

(define (random-unit)                   ; a flonum in [-1, 1)
  (set! seed (mod (+ (* 1103515245 seed) 12345) 2147483648))
  (- (/ (inexact seed) 1073741824.0) 1.0))

;; dense:   every coefficient random
;; zeros:   about a third of them 0.0 or -0.0, which the products skip
;; special: dense, but the scalar is +inf.0 and the last blade +nan.0
(define (coefficient kind)
  (let ((u (random-unit)))
    (cond ((and (eq? kind 'zeros) (< u -0.6)) 0.0)
          ((and (eq? kind 'zeros) (< u -0.3)) -0.0)
          (else u))))

(define (input p q r kind)
  (let* ((size (expt 2 (+ p q r)))
         (cs (let loop ((i 0) (acc '()))
               (if (= i size)
                   (reverse acc)
                   (loop (+ i 1) (cons (coefficient kind) acc))))))
    (list->multivector
      (make-algebra p q r)
      (if (eq? kind 'special)
          (append (list +inf.0) (cdr (reverse (cdr (reverse cs)))) (list +nan.0))
          cs))))

(define (emit op sig kind values)
  (write (append (list op sig kind) values))
  (newline))

(define (emit-mv op sig kind m)
  (emit op sig kind (multivector->list m)))

(for-each
  (lambda (sig)
    (for-each
      (lambda (kind)
        (let* ((p (car sig)) (q (cadr sig)) (r (caddr sig))
               (a (input p q r kind))
               (b (input p q r kind))
               (c (input p q r kind)))
          (emit-mv 'mv* sig kind (mv* a b))
          (emit-mv 'mv*3 sig kind (mv* a b c))
          (emit-mv 'mv-wedge sig kind (mv-wedge a b))
          (emit-mv 'mv-lcontract sig kind (mv-lcontract a b))
          (emit-mv 'mv+ sig kind (mv+ a b))
          (emit-mv 'mv- sig kind (mv- a b))
          (emit-mv 'mv-negate sig kind (mv- a))
          (emit-mv 'mv-scale sig kind (mv-scale a 0.3))
          (emit-mv 'mv-reverse sig kind (mv-reverse a))
          (emit-mv 'mv-involute sig kind (mv-involute a))
          (emit-mv 'mv-conjugate sig kind (mv-conjugate a))
          (emit-mv 'mv-dual sig kind (mv-dual a))
          (let loop ((k 0))
            (when (<= k (+ p q r))
              (emit-mv 'mv-grade sig kind (mv-grade a k))
              (loop (+ k 1))))
          (emit 'mv-norm2 sig kind (list (mv-norm2 a)))
          (emit 'mv-norm sig kind (list (mv-norm a)))
          ;; Not (= n2 0.0) alone: Curry's = is true when either side is
          ;; NaN.  Both mv-normalizes raise only on exactly zero.
          (unless (let ((n2 (mv-norm2 a))) (and (not (nan? n2)) (= n2 0.0)))
            (emit-mv 'mv-normalize sig kind (mv-normalize a)))))
      '(dense zeros special)))
  '((1 0 0) (0 1 0) (0 0 1) (2 0 0) (3 0 0) (0 3 0) (1 3 0)
    (3 0 1) (4 1 0) (2 2 2) (3 3 2) (8 0 0)))
