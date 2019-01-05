(in-package #:bayesian-analysis)

;; fixme: this needs to be refactored to be able to handle "inverse priors"

;; fixme: this needs some serious optimization, there is way too much
;; unoptimized stuff in here


(declaim (optimize (debug 3) (safety 1) (speed 3)))


(defgeneric prior-log-lambda (prior-type model slot-name)
  (:method ((type (eql :jeffreys)) model slot-name)
    (let+ ((min (get-value-of-slot/cat :min slot-name model))
	   (max (get-value-of-slot/cat :max slot-name model)))
      (if (or (<= min 0d0)
	      (<= max 0d0))
	  (error 'parameter-out-of-range
		 :parameter slot-name
		 :description
		 (format nil "Jeffrey priors are only applicable to (non-negative) scale parameters. 
The given prior range is: [~f,~f]" min max)))
      (let+ ((ln-max/min (coerce (log (/ max min)) 'double-float)))
	#'(lambda () (log (/ 1 (* (slot-value model slot-name) ln-max/min)))))))
  (:method ((type (eql :uniform)) model slot-name)
    (let+ ((min (get-value-of-slot/cat :min slot-name model))
	   (max (get-value-of-slot/cat :max slot-name model))
	   (1/max-min (coerce (/ 1 (- max min)) 'double-float)))
      #'(lambda () (log 1/max-min))))
  (:method ((type (eql :certain)) model slot-name)
    #'(lambda () 0d0)))



(defgeneric prior-is-constant (prior-type)
  (:method ((type (eql :certain))) t)
  (:method ((type (eql :uniform))) t)
  (:method ((type t)) nil))


(defclass gaussian-prior ()
  ((sigma :accessor sigma :initarg :sigma :initform 1d0)
   (mu :initarg :mu :accessor mu 
       :initform (error "Must initialize mu."))))


(defmethod prior-log-lambda ((prior gaussian-prior) model slot-name)
  (let+ (((&slots sigma mu) prior)
	 (2s^2 (* 2d0 sigma sigma))
	 (factor (log (/ 1d0 (* sigma (sqrt (* 2d0 pi)))))))
    #'(lambda ()
	(let* ((x-mu (- (slot-value model slot-name) mu))
	       (exponent (- (/ (* x-mu x-mu) 2s^2))))
	  (+ factor exponent)))))


;; ---

(defun make-prior-in-range (object slot-name)
  (let+ ((min (get-value-of-slot/cat :min slot-name object))
	 (max (get-value-of-slot/cat :max slot-name object)))
    #'(lambda ()
	(declare (type double-float min max))
	(<= min (slot-value object slot-name) max))))



;;;; model object initialization logic for priors

(defun --init--params-to-marginalize (model all-model-parameters)
  (iter
    (for name in all-model-parameters)
    (for marginalize-slot-name
	 = (w/suffix-slot-category :marginalize name))
    (if (slot-value model marginalize-slot-name)
	(collect name))))


(defun --init--prior-in-range-array (object model-parameters-to-marginalize)
  (let+ ((no-priors (length model-parameters-to-marginalize))
	 (prior-in-range-array (make-array no-priors :element-type '(function () t)
						     :initial-element #'(lambda () nil))))
    (iter
      (for slot-name in-sequence model-parameters-to-marginalize with-index i)
      (for in-range = (make-prior-in-range  object slot-name))
      (setf (aref prior-in-range-array i) in-range)
      (finally (return prior-in-range-array)))))


(defun --init--prior-array (object model-params-to-marginalize)
  (let+ ((varying-priors-slots
	  (remove-if #'prior-is-constant model-params-to-marginalize
		     :key #'(lambda (slot-name)
			      (get-value-of-slot/cat :prior slot-name object))))
	 (varying/log-priors-array (make-array (length varying-priors-slots)
					       :element-type '(function () t)
					       :initial-element #'(lambda () 0d0))))
    (iter
      (for slot-name in-sequence varying-priors-slots with-index i)
      (for type = (get-value-of-slot/cat :prior slot-name object))
      (for prior = (prior-log-lambda type object slot-name))
      (setf (aref varying/log-priors-array i) prior)
      (finally (return (values varying/log-priors-array (length varying/log-priors-array)))))))


(defun --init--all-priors-array (object parameters)
  (let+ ((varying/log-priors-array (make-array (length parameters)
					       :element-type '(function () t)
					       :initial-element #'(lambda () 0d0))))
    (iter
      (for slot-name in-sequence parameters with-index i)
      (for type = (get-value-of-slot/cat :prior slot-name object))
      (setf (aref varying/log-priors-array i)
	    (prior-log-lambda type object slot-name))
      (finally (return (values varying/log-priors-array (length varying/log-priors-array)))))))


(defun --init--get-log-constant-prior (object model-params-to-marginalize)
  (let+ ((constant-priors-slots
	  (remove-if-not #'prior-is-constant model-params-to-marginalize
			 :key #'(lambda (slot-name)
				  (get-value-of-slot/cat :prior slot-name object))))
	 (log-of-constant-priors (reduce #'+ constant-priors-slots
					 :key #'(lambda (slot-name)
						  (funcall
						   (prior-log-lambda (get-value-of-slot/cat :prior slot-name object)
								     object
								     slot-name))))))
    log-of-constant-priors))


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







