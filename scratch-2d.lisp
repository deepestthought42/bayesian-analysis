(in-package #:bayesian-analysis)


(declaim (optimize (debug 3) (space 0) (safety 1) (speed 3)))



(ba:define-data-class 1d-gaussian x y 0
    (obj (source t))
  (let+ ((no-data-points 2000)
	 (arr (gsl-cffi:get-array-random-uniform no-data-points))
	 (mean1 -2d0)
	 (mean2 3d0)
	 (a1 10)
	 (a2 1)
	 (sigma 0.45d0)
	 (rng (gsl-cffi:get-random-number-generator gsl-cffi::*mt19937_1999* 11)))
    #+sbcl (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
    (setf x (make-array no-data-points
			:element-type 'double-float
			:initial-contents
			(iter
			  (for i from 0 below no-data-points)
			  (collect (if (< (/ a2 a1) (aref arr i))
				       (+ mean1 (gsl-cffi:random-gaussian rng sigma))
				       (+ mean2 (gsl-cffi:random-gaussian rng sigma))))))
	  y (make-array no-data-points
			:element-type 'double-float
			:initial-element 1d0)
	  ;; err (make-array no-data-points
	  ;; 		  :element-type 'double-float
	  ;; 		  :initial-element 0.1d0)
	  )))


(define-bayesian-model (test-mean 1d-gaussian)
    ((mu-1 :marginalize t :prior-type :uniform :min -40 :max 40 :default -2d0 :sample-sigma 0.01d0)
     (mu-2 :marginalize t :prior-type :uniform :min -40 :max 40 :default 3d0 :sample-sigma 0.01d0)
     (A-1 :marginalize t :prior-type :jeffreys :min 0.01 :max 2 :default 1 :sample-sigma 0.001d0)
     (A-2 :marginalize t :prior-type :jeffreys :min 0.01 :max 2 :default 0.1 :sample-sigma 0.001d0)
     (sigma :marginalize t :prior-type :jeffreys :min 0.01 :max 10 :default 0.5 :sample-sigma 0.01d0))
    (:p_of_x_i=f_i_of_x_i)
    ((x)
      (let* ((A (+ A-1 A-2))
	     (e1 (*  (/ A-1 A)
		     (exp (- (/ (expt (- x mu-1) 2d0)
				( * 2d0 sigma sigma))))))
	     (e2 (*  (/ A-2 A)
		     (exp (- (/ (expt (- x mu-2) 2d0)
				( * 2d0 sigma sigma))))))
	     (val (/ (+ e1 e2)
		     (* (sqrt (* 2d0 pi)) sigma))))
	(declare (type (double-float 0d0) e1 e2 a val))
	(if (= 0d0 val)
	    most-negative-double-float
	     val))))



