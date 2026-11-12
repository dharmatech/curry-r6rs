#!r6rs
;;; Compare two files of data written by the differential tests.
;;;
;;;   scheme --program tests/compare-results.sps expected.txt actual.txt
;;;
;;; Numbers must be eqv?, so 0.0 and -0.0 differ.  Any NaN matches any NaN:
;;; IEEE leaves NaN sign and payload unspecified, and x86 picks them by
;;; operand order.  Exits 0 when everything matches.

(import (rnrs))

(define (same? x y)
  (cond ((and (number? x) (number? y))
         (or (eqv? x y)
             (and (real? x) (real? y) (nan? x) (nan? y))))
        ((and (pair? x) (pair? y))
         (and (same? (car x) (car y)) (same? (cdr x) (cdr y))))
        (else (equal? x y))))

(define (read-all file)
  (call-with-input-file file
    (lambda (port)
      (let loop ((acc '()))
        (let ((x (read port)))
          (if (eof-object? x) (reverse acc) (loop (cons x acc))))))))

(define args (cdr (command-line)))
(unless (= (length args) 2)
  (display "usage: compare-results.sps expected actual\n" (current-error-port))
  (exit 2))

(define expected (read-all (car args)))
(define actual (read-all (cadr args)))

(define (count-numbers x)
  (cond ((number? x) 1)
        ((pair? x) (+ (count-numbers (car x)) (count-numbers (cdr x))))
        (else 0)))

(define mismatches
  (let loop ((e expected) (a actual) (n 0))
    (cond ((or (null? e) (null? a)) n)
          ((same? (car e) (car a)) (loop (cdr e) (cdr a) n))
          (else
           (when (< n 5)
             (display "MISMATCH\n  expected: ") (write (car e))
             (display "\n  actual:   ") (write (car a)) (newline))
           (loop (cdr e) (cdr a) (+ n 1))))))

(unless (= (length expected) (length actual))
  (display "different number of results: ")
  (display (length expected)) (display " vs ") (display (length actual)) (newline))

(display (length expected)) (display " results, ")
(display (count-numbers expected)) (display " numbers, ")
(display mismatches) (display " mismatched") (newline)

(exit (if (and (= mismatches 0) (= (length expected) (length actual))) 0 1))
