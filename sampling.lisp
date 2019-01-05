(in-package #:bayesian-analysis)


(defgeneric get-sampler (sampler-type ))



(defclass sampler ()
  ((rng :initarg :rng :accessor rng 
	:initform (gsl-cffi:get-random-number-generator gsl-cffi::*mt19937_1999*))))

(defclass gaussian-sampler (sampler)
  ((sigma :accessor sigma :initarg :sigma :initform 0.1d0)))

(defgeneric get-sampling-lambda (sampler object slot-name))

(defmethod get-sampling-lambda ((sampler gaussian-sampler) object slot-name)
  (let+ (((&slots rng sigma) sampler))
    #'(lambda ()
	(setf (slot-value object slot-name)
	      (+ (slot-value object slot-name)
		 (gsl-cffi:random-gaussian rng sigma))))))



;;; model object initialization code goes here


(defun --init--sampling-fn-array (object parameters-to-marginalize)
  (let+ ((no-samplers (1+ (length parameters-to-marginalize)))
	 (sampler-array (make-array no-samplers :element-type '(function ())
						:initial-element #'(lambda () (setf (new-sample? object) t)))))
    (iter
      (for slot-name in-sequence parameters-to-marginalize with-index i)
      (for slot-val = (get-value-of-slot/cat :sampler slot-name object))
      (for slot-sampler-sigma = (get-value-of-slot/cat :sample-sigma slot-name object))
      (for sampler = (if (symbolp slot-val) (make-instance slot-val :sigma slot-sampler-sigma) slot-val))
      (for sampling-lambda = (get-sampling-lambda sampler object slot-name))
      (setf (aref sampler-array i) sampling-lambda)
      (finally (return (values sampler-array no-samplers))))))



(defun --init--sampling-functions (model-object)
  (let+ (((&slots model-parameters-to-marginalize
		  sample-new-parameters) model-object)
	 ((&values array-sampling-fun no-samplers)
	  (--init--sampling-fn-array model-object model-parameters-to-marginalize)))
    (labels ((sample-new-params ()
	       (iter
    		 (for i from 0 below no-samplers)
    		 (funcall (aref array-sampling-fun i))
    		 (finally (return t)))))
      (setf sample-new-parameters #'sample-new-params))))



