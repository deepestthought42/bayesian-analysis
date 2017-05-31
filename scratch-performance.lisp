(in-package #:fit-penning)


(declaim (optimize (debug 0) (safety 0) (speed 3) (space 0)))

(defmethod create-parameter-array ((model (eql 'fp:2d-bayes-konig)))
  (static-vectors:make-static-vector 19 :element-type 'double-float
				     :initial-element 0d0))


(defstruct 2d-bayes-konig-parameters
  ;; dependent variables
  (f-of-x 0d0 :type double-float)
  ;; cache
  (nu-rf-cache 0d0 :type (double-float 0d0))
  (tof-cache 0d0 :type (double-float 0d0))
  (f-of-x-cache 0d0 :type double-float)
  (new-sample t :type boolean))


(defmethod array-to-values ((model fp:2d-bayes-konig) array)
  (declare (type (simple-array double-float (19)) array))
  (setf (om-c model) (aref array 0)
	(q model) (aref array 1)
	(om-m model) (aref array 2)
	(om-rf model) (aref array 3)
	(b-max model) (aref array 4)
	(e-0 model) (aref array 5)
	(v-0 model) (aref array 6)
	(t-rf model) (aref array 7)
	(no-conv model) (aref array 8)
	(rho-m0 model) (aref array 9)
	(rho-p0 model) (aref array 10)
	(damping model) (aref array 11)
	(delta-phi model) (aref array 12)
	(tof-offset model) (aref array 13)
	(sigma model) (aref array 14)
	(p-interest model) (aref array 15)
	(contaminant-offset model) (aref array 16)
	(z-start model) (aref array 17)
	(z-end model) (aref array 18))
  model)


(defmethod values-to-array ((model fp:2d-bayes-konig) &key (array (create-parameter-array 'fp:2d-bayes-konig)))
  (declare (type (simple-array double-float (19)) array))
  (let+ (((&slots om-c q om-m om-rf b-max e-0 v-0 t-rf no-conv rho-m0 rho-p0
		  damping delta-phi tof-offset sigma p-interest
		  contaminant-offset z-start z-end) model))
    (setf (aref array 0) (om-c model)
	  (aref array 1) (q model)
	  (aref array 2) (om-m model)
	  (aref array 3) (om-rf model)
	  (aref array 4) (b-max model)
	  (aref array 5) (e-0 model)
	  (aref array 6) (v-0 model)
	  (aref array 7) (t-rf model)
	  (aref array 8) (no-conv model)
	  (aref array 9) (rho-m0 model)
	  (aref array 10) (damping model)
	  (aref array 11) (rho-p0 model)
	  (aref array 12) (delta-phi model)
	  (aref array 13) (tof-offset model)
	  (aref array 14) (sigma model)
	  (aref array 15) (p-interest model)
	  (aref array 16) (contaminant-offset model)
	  (aref array 17) (z-start model)
	  (aref array 18) (z-end model))
    array))

(declaim (inline mass)
	 (ftype (function (double-float double-float double-float) double-float) mass))
(defun mass (om-c b-max q)
  (declare (type double-float om-c b-max q))
  (the (double-float 0.0d0) (/ (* q b-max) om-c)))


(declaim (inline integrand)
	 (ftype (function (double-float double-float double-float double-float double-float
					  double-float double-float double-float double-float
					  double-float double-float double-float double-float
					  double-float)
			    double-float) integrand))
(defun integrand (om-rf om-m om-c b-max delta-phi q e-0 v-0
		  no-conv damping rho-m0 rho-p0 t-rf z-start)
  (declare (type double-float om-rf om-m om-c b-max delta-phi q e-0 v-0
		 no-conv damping rho-m0 rho-p0 t-rf)
	   (ftype (function (double-float double-float double-float) double-float)
		  interpolated-integration)
	   (inline interpolated-integration)
	   (ignorable delta-phi v-0 damping rho-p0))
  (let ((e0/q (/ e-0 q))
	(mu/q
	  (/
	   (%mu-over-bmax nil nil om-rf om-c om-m rho-m0 rho-p0
			  damping no-conv t-rf q :have-delta-phi
			  nil)
	   q)))
    (* (sqrt (the (double-float 0.0d0) (/ (mass om-c b-max q) q)))
       (interpolated-integration e0/q mu/q z-start))))


(declaim (inline mean)
	 (ftype (function (double-float double-float double-float double-float double-float
					double-float double-float double-float double-float
					double-float double-float double-float double-float
					double-float double-float)
			  double-float)
		mean))

(defun mean (om-rf tof-offset om-m om-c b-max delta-phi q e-0 v-0
	     no-conv damping rho-m0 rho-p0 t-rf z-start)
  (declare (type (function (double-float) double-float) e-spline b-spline)
	   (type double-float om-rf tof-offset om-m om-c b-max delta-phi q e-0 v-0
		 no-conv damping rho-m0 rho-p0 t-rf z-start))
  (+ tof-offset (* 1000000.0d0 (integrand om-rf om-m om-c b-max delta-phi q e-0 v-0
					  no-conv damping rho-m0 rho-p0 t-rf z-start))))

(declaim (inline standard-normal)
	 (ftype (function (double-float) double-float)))

(defun standard-normal (x)
  (declare (type double-float x))
  (the double-float
       (/ (exp (/ (- (* x x)) 2.0d0)) (sqrt (* 2.0d0 pi)))))

(declaim (inline gaussian)
	 (ftype (function (double-float double-float double-float double-float) double-float)
		gaussian))
(defun gaussian (intensity x x-0 sigma)
  (declare (type double-float x x-0 sigma intensity))
  (the double-float
       (* intensity
	  (/ (standard-normal (/ (- x x-0) sigma)) sigma))))

(declaim (inline cached-mean)
	 (ftype (function (double-float double-float double-float double-float 
					double-float double-float double-float double-float
					double-float double-float double-float double-float
					double-float double-float double-float 2d-bayes-konig-parameters)
			  double-float) cached-mean))

(defun cached-mean (nu-rf om-m om-c b-max delta-phi q e-0 v-0
		    no-conv damping rho-m0 rho-p0 t-rf z-start tof-offset parameters)
  (declare (type 2d-bayes-konig-parameters parameters)
	   (type double-float nu-rf tof-offset om-m om-c b-max delta-phi q e-0 v-0
		 no-conv damping rho-m0 rho-p0 t-rf z-start))
  (let+ (((&structure 2d-bayes-konig-parameters- nu-rf-cache
		      f-of-x-cache new-sample) parameters))
    (if (and (not new-sample)
	     (= nu-rf-cache nu-rf))
	f-of-x-cache
	(setf new-sample nil
	      nu-rf-cache nu-rf
	      f-of-x-cache
	      (mean (* 2.0d0 pi nu-rf) tof-offset om-m om-c b-max delta-phi q e-0 v-0
		    no-conv damping rho-m0 rho-p0 t-rf z-start)))))



(defun perf-test/2d-bayes-konig-model-function (nu-rf tof parameters array)
  (declare (type 2d-bayes-konig-parameters parameters)
	   (type (simple-array double-float (19)) array)
	   (type (double-float 0d0) nu-rf tof))

  (let+ ((om-c (aref array 0))
	 (q (aref array 1))
	 (om-m (aref array 2))
	 ;(om-rf (aref array 3))
	 (b-max (aref array 4))
	 (e-0 (aref array 5))
	 (v-0 (aref array 6))
	 (t-rf (aref array 7))
	 (no-conv (aref array 8))
	 (rho-m0 (aref array 9))
	 (rho-p0 (aref array 10))
	 (damping (aref array 11))
	 (delta-phi (aref array 12))
	 (tof-offset (aref array 13))
	 (sigma (aref array 14))
	 (p-interest (aref array 15))
	 (contaminant-offset (aref array 16))
	 (z-start (aref array 17))
	 ;(z-end (aref array 18))
	 ((&structure 2d-bayes-konig-parameters- f-of-x) parameters))
    (declare (type double-float om-c q om-m  b-max e-0 v-0 t-rf no-conv rho-m0
		   rho-p0 damping delta-phi tof-offset sigma p-interest contaminant-offset
		   z-start))
    (declare (type (double-float 0.0d0) tof nu-rf))
    (setf f-of-x
	  (+ (gaussian p-interest
		       tof
		       (cached-mean nu-rf om-m om-c b-max delta-phi q e-0 v-0
				    no-conv damping rho-m0 rho-p0 t-rf z-start
				    tof-offset
				    parameters)
		       sigma)
	     (gaussian (- 1.0d0 p-interest) tof contaminant-offset sigma)))))





(defun test-function ()
  (let+ ((in/out-val (make-2d-bayes-konig-parameters))
	 (model (pa::f/make-model))
	 (array (values-to-array model)))
    (time
     (dotimes (v 10000000)
       (perf-test/2d-bayes-konig-model-function 1000000d0 30d0 in/out-val array)))))



(defun test-function2 ()
  (let+ ((in/out-val (make-2d-bayes-konig-parameters))
	 (model (pa::f/make-model))
	 (array (values-to-array model)))
    (time
     (dotimes (v 10000000)
       (2d-bayes-konig-model-function 1000000d0 30d0 model)))))




