#!r6rs
;;; Benchmark: (curry-r6rs multivector) on Chez Scheme.
;;;
;;; bench/multivector-curry-c.scm runs the same workloads against Curry's
;;; built-in C multivectors.  Keep the two files' workloads in step: same
;;; inputs, same order of operations, so both print the same checksums.
;;;
;;; Only the timer and the collector call come from (chezscheme); everything
;;; that is measured is portable R6RS.
;;;
;;; Run from the repository root:  scheme --program bench/multivector.sps

(import (rnrs)
        (curry-r6rs multivector)
        (only (chezscheme) current-time time-second time-nanosecond collect))

(define repetitions 5)

(define (now-ns)
  (let ((t (current-time 'time-monotonic)))
    (+ (* (time-second t) 1000000000) (time-nanosecond t))))

;;; ---------------------------------------------------------------------------
;;; Workloads -- keep in step with multivector-curry-c.scm

;; A multivector with every coefficient nonzero, so no product term is skipped.
(define (dense-mv alg seed)
  (let loop ((i (- (algebra-blade-count alg) 1)) (acc '()))
    (if (< i 0)
        (list->multivector alg acc)
        (loop (- i 1)
              (cons (let ((c (/ (inexact (+ 1 (mod (+ (* i 7) seed) 11))) 7.0)))
                      (if (odd? i) (- c) c))
                    acc)))))

;; The geometric product of two dense multivectors, over and over.
(define (product-workload p q r iterations)
  (let* ((alg (make-algebra p q r))
         (a (dense-mv alg 1))
         (b (dense-mv alg 2)))
    (lambda ()
      (let loop ((k 0) (sum 0.0))
        (if (= k iterations)
            sum
            (loop (+ k 1) (+ sum (mv-scalar (mv* a b)))))))))

;; Rotate `count` vectors by a fixed rotor, R v R~, building each vector
;; from a list the way application code would.
(define (rotor-workload count)
  (let* ((g3 (make-algebra 3 0 0))
         (len (sqrt 14.0))
         (ax (/ 1.0 len)) (ay (/ 2.0 len)) (az (/ 3.0 len))
         (h 0.35)
         (s (- (sin h)))
         (r (list->multivector g3 (list (cos h) 0.0 0.0 (* s az) 0.0 (* (- s) ay) (* s ax) 0.0)))
         (rrev (mv-reverse r)))
    (lambda ()
      (let loop ((i 0) (sum 0.0))
        (if (= i count)
            sum
            (let* ((x (* (inexact (mod i 1000)) 0.001))
                   (v (list->multivector g3 (list 0.0 x (- 1.0 x) 0.0 0.5 0.0 0.0 0.0)))
                   (w (mv* (mv* r v) rrev)))
              (loop (+ i 1)
                    (+ (+ (+ sum (multivector-ref w 1))
                          (multivector-ref w 2))
                       (multivector-ref w 4)))))))))

;;; ---------------------------------------------------------------------------
;;; Harness

(define (insert x sorted)
  (cond ((null? sorted) (list x))
        ((<= x (car sorted)) (cons x sorted))
        (else (cons (car sorted) (insert x (cdr sorted))))))

(define (ms ns) (/ (inexact ns) 1e6))

;; Prints: implementation workload algebra operations best-ms median-ms
;;         ns-per-operation checksum
(define (bench name algebra operations thunk)
  (let loop ((k 0) (times '()) (checksum #f))
    (if (= k repetitions)
        (let ((best (car times))
              (median (list-ref times (div repetitions 2))))
          (for-each (lambda (x) (display x) (display " "))
                    (list "chez" name algebra operations (ms best) (ms median)
                          (/ (inexact best) operations)))
          (display checksum)
          (newline))
        (begin
          (collect)
          (let* ((t0 (now-ns))
                 (result (thunk))
                 (t1 (now-ns)))
            (loop (+ k 1) (insert (- t1 t0) times) result))))))

(bench "product" "Cl(3,0,0)" 1000000 (product-workload 3 0 0 1000000))
(bench "product" "Cl(4,1,0)" 200000  (product-workload 4 1 0 200000))
(bench "product" "Cl(8,0,0)" 1000    (product-workload 8 0 0 1000))
(bench "rotor"   "Cl(3,0,0)" 500000  (rotor-workload 500000))
