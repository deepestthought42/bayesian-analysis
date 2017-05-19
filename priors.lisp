(in-package #:bayesian-analysis)

;; fixme: this needs to be refactored to be able to handle "inverse priors"

;; fixme: this needs some serious optimization, there is way too much
;; unoptimized stuff in here


(declaim (optimize (debug 3) (safety 1) (speed 3)))



;; --- jeffreys prior

(defun jeffreys-log-lambda (min max object slot-name)
  (let+ ((ln-max/min (coerce (log (/ max min)) 'double-float)))
    (if *debug-function*
	#'(lambda ()
	    (let ((val (log (/ 1 (* (slot-value object slot-name) ln-max/min)))))
	      (debug-out :info :prior
			 "Getting Jeffrey prior value for slot: ~a, value: ~f"
			 slot-name val)
	      val))
	#'(lambda ()
	    (log (/ 1 (* (slot-value object slot-name) ln-max/min)))))))



(defun jeffreys (min max object slot-name)
  (let ((ln-max/min (coerce (log (/ max min)) 'double-float)))
    (declare (type double-float ln-max/min))
    (log (/ 1 (* (slot-value object slot-name) ln-max/min)))))


;; --- uniform prior

(defun uniform-log-lambda (min max object slot-name)
  (declare (ignore object))
  (let ((1/delta-x (coerce (/ 1 (- max min)) 'double-float)))
    (if *debug-function*
	#'(lambda ()
	    (declare (type double-float 1/delta-x))
	    (debug-out :info :prior
		       "Getting uniform prior value for slot: ~a, value: ~f"
		       slot-name 1/delta-x)
	    (log 1/delta-x))
	#'(lambda ()
	    (declare (type double-float 1/delta-x))
	    (log 1/delta-x)))))

(defun uniform (min max object slot-name)
  (declare (ignore object slot-name))
  (coerce (/ 1 (- max min)) 'double-float))

;; --- certain prior

(defun certain-log-lambda (min max object slot-name)
  (declare (ignore object min max))
  (if *debug-function*
      #'(lambda ()
	  (debug-out :info :prior
		     "Getting certain prior value for slot: ~a"
		     slot-name)
	  0d0)
      #'(lambda () 0d0)))

(defun certain (min max object slot-name)
  (declare (ignore object slot-name min max))
  1d0)

;; --- gaussian prior

(defun gaussian-log-lambda (min max object slot-name)
  (declare (ignore min max))
  (labels ((sv (cat)
	     (slot-value object (w/suffix-slot-category cat slot-name))))
    (let+ ((sigma (sv :prior-sigma))
	   (2s^2 (* 2d0 sigma sigma))
	   (mu (sv :prior-mu))
	   (factor (log (/ 1d0 (* sigma (sqrt (* 2d0 pi)))))))
      (if *debug-function*
	  #'(lambda ()
	      (let* ((x-mu (- (slot-value object slot-name) mu))
		     (exponent (/ (* x-mu x-mu) 2s^2))
		     (val (+ factor exponent)))
		(debug-out :info :prior
			   "Getting Gaussian prior value for slot: ~a, x: ~a, value: ~f"
			   slot-name (slot-value object slot-name) val)
		val))
	  #'(lambda ()
	      (let* ((x-mu (- (slot-value object slot-name) mu))
		     (exponent (/ (* x-mu x-mu) 2s^2)))
		(+ factor exponent)))))))


(defun gaussian (min max object slot-name)
  (declare (ignore min max))
  (labels ((sv (cat)
	     (slot-value object (w/suffix-slot-category cat slot-name))))
    (let+ ((sigma (sv :prior-sigma))
	   (2s^2 (* 2d0 sigma sigma))
	   (mu (sv :prior-mu))
	   (factor (log (/ 1d0 (* sigma (sqrt (* 2d0 pi))))))
	   (x-mu (- (slot-value object slot-name) mu))
	   (exponent (/ (* x-mu x-mu) 2s^2)))
      (+ factor exponent))))


;; ---

(defun make-prior-in-range (min max object slot-name)
  (if *debug-function*
      #'(lambda ()
	  (declare (type double-float min max))
	  (let+ ((val (slot-value object slot-name))
		 (in-range (<= min val max)))
	    (debug-out :info :prior
		       "Is prior for slot: ~a and value: ~f in range: ~a"
		       slot-name val in-range)
	    in-range))
      #'(lambda ()
	  (declare (type double-float min max))
	  (<= min (slot-value object slot-name) max))))


(defparameter *prior-types*
  '((:jeffreys jeffreys-log-lambda :variable jeffreys)
    (:uniform uniform-log-lambda :constant uniform)
    (:certain certain-log-lambda :constant certain)
    (:gaussian gaussian-log-lambda :variable)))



(defun find-prior-type (type &key (type-list *prior-types*))
  (alexandria:if-let (type (find type type-list :key #'first))
    type (error 'unknown-prior-type :type-not-known type)))

(defun prior-is-constant (type &key (type-list *prior-types*))
  (equal :constant (third (find-prior-type type :type-list type-list))))

(defun get-prior-lambda (type &key (type-list *prior-types*))
  (second (find-prior-type type :type-list type-list)))

(defun get-prior-function (type &key (type-list *prior-types*))
  (fourth (find-prior-type type :type-list type-list)))

(defun get-inverse-prior-lambda (type &key (type-list *prior-types*))
  (lambda () ))



;;;; model object initialization logic for priors

(defun --init--params-to-marginalize (model all-model-parameters)
  (iter
    (for name in all-model-parameters)
    (for marginalize-slot-name
	 = (w/suffix-slot-category :marginalize name))
    (if (slot-value model marginalize-slot-name)
	(collect name))))


(defun --init--prior-in-range-array (object model-parameters-to-marginalize)
  (labels ((sn (cat name)
	     (w/suffix-slot-category cat name))
	   (sv (cat name)
	     (slot-value object (sn cat name))))
    (let+ ((no-priors (length model-parameters-to-marginalize))
	   (prior-in-range-array (make-array no-priors :element-type '(function () t)
						       :initial-element #'(lambda () nil))))
      (iter
	(for slot-name in-sequence model-parameters-to-marginalize with-index i)
	(for in-range = (make-prior-in-range (sv :min slot-name) (sv :max slot-name) object slot-name))
	(setf (aref prior-in-range-array i) in-range)
	(finally (return prior-in-range-array))))))


(defun --init--prior-array (object model-params-to-marginalize)
  (labels ((sn (cat name)
	     (w/suffix-slot-category cat name))
	   (sv (cat name)
	     (slot-value object (sn cat name))))
    (let+ ((varying-priors-slots
	    (remove-if #'prior-is-constant model-params-to-marginalize
		       :key #'(lambda (slot-name) (sv :prior-type slot-name))))
	   (varying/log-priors-array (make-array (length varying-priors-slots)
						 :element-type '(function () t)
						 :initial-element #'(lambda () 0d0))))
      (iter
	(for slot-name in-sequence varying-priors-slots with-index i)
	(for prior = (funcall (get-prior-lambda (sv :prior-type slot-name))
			      (sv :min slot-name)
			      (sv :max slot-name)
			      object slot-name))
	(setf (aref varying/log-priors-array i) prior)
	(finally (return (values varying/log-priors-array (length varying/log-priors-array))))))))


(defun --init--all-priors-array (object parameters)
  (labels ((sn (cat name)
	     (w/suffix-slot-category cat name))
	   (sv (cat name)
	     (slot-value object (sn cat name))))
    (let+ ((varying/log-priors-array (make-array (length parameters)
						 :element-type '(function () t)
						 :initial-element #'(lambda () 0d0))))
      (iter
	(for slot-name in-sequence parameters with-index i)
	(setf (aref varying/log-priors-array i)
	      (funcall (get-prior-lambda (sv :prior-type slot-name))
		       (sv :min slot-name)
		       (sv :max slot-name)
		       object slot-name))
	(finally (return (values varying/log-priors-array (length varying/log-priors-array))))))))


(defun --init--get-log-constant-prior (object model-params-to-marginalize)
    (labels ((sn (cat name)
	     (w/suffix-slot-category cat name))
	   (sv (cat name)
	     (slot-value object (sn cat name))))
    (let+ ((constant-priors-slots
	    (remove-if-not #'prior-is-constant model-params-to-marginalize
			   :key #'(lambda (slot-name)
				    (slot-value object (sn :prior-type slot-name)))))
	   (log-of-constant-priors (reduce #'+ constant-priors-slots
					   :key #'(lambda (slot-name)
						    (log (funcall (get-prior-function (sv :prior-type slot-name))
								  (sv :min slot-name) (sv :max slot-name)
								  object
								  slot-name))))))

      log-of-constant-priors)))


(defun --init--priors (model-object)
  (let+ (((&slots priors-in-range constant/log-of-priors varying/log-of-priors
		  model-parameters-to-marginalize model-prior all-model-parameters
		  log-of-all-priors log-of-all-priors-array)
	  model-object)
	 (no-priors (length model-parameters-to-marginalize))
	 (log-of-constant-priors (--init--get-log-constant-prior model-object model-parameters-to-marginalize))
	 ((&values array-of-log-of-varying-priros no-varying-priors)
	  (--init--prior-array model-object model-parameters-to-marginalize))
	 ((&values array-of-log-of-all-priors no-of-all-priors)
	  (--init--all-priors-array model-object all-model-parameters))
	 (array-of-prior-in-range (--init--prior-in-range-array model-object model-parameters-to-marginalize)))
    (labels ((priors-in-range ()
    	       (iter
    		 (for i from 0 below no-priors)
    		 (if (not (funcall (aref array-of-prior-in-range i)))
    		     (return nil))
    		 (finally (return t))))
    	     (log-of-varying-prior ()
    	       (iter
		 (declare (type double-float log-of-p))
    		 (for i from 0 below no-varying-priors)
    		 (iter:multiply (funcall (aref array-of-log-of-varying-priros i))
				into log-of-p)
    		 (finally (return log-of-p))))
	     ;; fixme: there might be room for optimization here: for
	     ;; example, one might consider to create a special
	     ;; function at compile time (macroexpansion time) for
	     ;; this
	     (log-of-all-priors ()
    	       (iter
    		 (for i from 0 below no-of-all-priors)
    		 (iter:multiply (funcall (aref array-of-log-of-all-priors i))
				into log-of-p)
    		 (finally (return log-of-p)))))
      (setf priors-in-range #'priors-in-range
	    log-of-all-priors #'log-of-all-priors
	    constant/log-of-priors (constantly log-of-constant-priors)
	    varying/log-of-priors #'log-of-varying-prior
	    log-of-all-priors-array array-of-log-of-all-priors))))







