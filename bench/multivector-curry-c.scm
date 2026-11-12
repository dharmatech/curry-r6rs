;;; Benchmark: Curry Scheme's built-in C multivectors (src/multivec.c).
;;;
;;; The same workloads as bench/multivector.sps, written against Curry's own
;;; API (make-mv, mv-from-list, mv-ref, ...).  Keep the two files in step:
;;; same inputs, same order of operations, so both print the same checksums.
;;;
;;; This file is Curry-only, not portable.  Run it with:
;;;   curry --clear-cache bench/multivector-curry-c.scm

(define repetitions 5)

(define (now-ns) (current-jiffy))       ; nanoseconds in Curry

;;; ---------------------------------------------------------------------------
;;; Workloads -- keep in step with multivector.sps

(define (dense-mv p q r seed)
  (let loop ((i (- (expt 2 (+ p q r)) 1)) (acc '()))
    (if (< i 0)
        (mv-from-list p q r acc)
        (loop (- i 1)
              (cons (let ((c (/ (exact->inexact (+ 1 (modulo (+ (* i 7) seed) 11))) 8.0)))
                      (if (odd? i) (- c) c))
                    acc)))))

(define (product-workload p q r iterations)
  (let ((a (dense-mv p q r 1))
        (b (dense-mv p q r 2)))
    (lambda ()
      (let loop ((k 0) (sum 0.0))
        (if (= k iterations)
            sum
            (loop (+ k 1) (+ sum (mv-scalar (mv* a b)))))))))

(define (rotor-workload count)
  (let* ((len (sqrt 14.0))
         (ax (/ 1.0 len)) (ay (/ 2.0 len)) (az (/ 3.0 len))
         (h 0.35)
         (s (- (sin h)))
         (r (mv-from-list 3 0 0 (list (cos h) 0.0 0.0 (* s az) 0.0 (* s ay) (* s ax) 0.0)))
         (rrev (mv-reverse r)))
    (lambda ()
      (let loop ((i 0) (sum 0.0))
        (if (= i count)
            sum
            (let* ((x (* (exact->inexact (modulo i 1000)) 0.001))
                   (v (mv-from-list 3 0 0 (list 0.0 x (- 1.0 x) 0.0 0.5 0.0 0.0 0.0)))
                   (w (mv* (mv* r v) rrev)))
              (loop (+ i 1)
                    (+ (+ (+ sum (mv-ref w 1))
                          (mv-ref w 2))
                       (mv-ref w 4)))))))))

;;; ---------------------------------------------------------------------------
;;; Harness

(define (insert x sorted)
  (cond ((null? sorted) (list x))
        ((<= x (car sorted)) (cons x sorted))
        (else (cons (car sorted) (insert x (cdr sorted))))))

(define (ms ns) (/ (exact->inexact ns) 1e6))

(define (bench name algebra operations thunk)
  (let loop ((k 0) (times '()) (checksum #f))
    (if (= k repetitions)
        (let ((best (car times))
              (median (list-ref times (quotient repetitions 2))))
          (for-each (lambda (x) (display x) (display " "))
                    (list "curry-c" name algebra operations (ms best) (ms median)
                          (/ (exact->inexact best) operations)))
          (display checksum)
          (newline))
        (begin
          (gc)
          (let* ((t0 (now-ns))
                 (result (thunk))
                 (t1 (now-ns)))
            (loop (+ k 1) (insert (- t1 t0) times) result))))))

(bench "product" "Cl(3,0,0)" 1000000 (product-workload 3 0 0 1000000))
(bench "product" "Cl(4,1,0)" 200000  (product-workload 4 1 0 200000))
(bench "product" "Cl(8,0,0)" 1000    (product-workload 8 0 0 1000))
(bench "rotor"   "Cl(3,0,0)" 500000  (rotor-workload 500000))
