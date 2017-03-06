(in-package #:bayesian-analysis)





(defun gaussian-lambda (rng sigma object slot-name)
  (if *debug-function*
      #'(lambda ()
	  (let+ ((old (slot-value object slot-name))
		 (new (+ old (gsl-cffi:random-gaussian rng sigma))))
	    (debug-out :info :sample "Gaussian draw for: ~a, old: ~,10d, new: ~,10d"
		       slot-name old new)
	    (setf (slot-value object slot-name) new)))
      #'(lambda ()
	  (setf (slot-value object slot-name)
		(+ (slot-value object slot-name)
		   (gsl-cffi:random-gaussian rng sigma))))))



(defparameter *sampling-types*
  '((:gaussian gaussian-lambda)))


(define-condition unknown-sampling-type (error)
  ((type-not-known :accessor type-not-known :initarg :type-not-known
		   :initform :unspecified)))


(defun find-sampling-type (type &key (type-list *sampling-types*))
  (alexandria:if-let (type (find type type-list :key #'first))
    type (error 'unknown-sampling-type :type-not-known type)))



(defun get-sampling-creator (type &key (type-list *sampling-types*))
  (second (find-sampling-type type :type-list type-list)))


;;; model object initialization code goes here


(defun --init--sampling-fn-array (object parameters-to-marginalize)
  (labels ((sn (cat name)
	     (w/suffix-slot-category cat name))
	   (sv (cat name)
	     (slot-value object (sn cat name))))
    (let+ ((no-samplers (length parameters-to-marginalize))
	   (sampler-array (make-array no-samplers :element-type '(function ())
						  :initial-element #'(lambda ())))
	   ((&slots rng) object))
      (iter
	(for slot-name in-sequence parameters-to-marginalize with-index i)
	(for sampler = (funcall (get-sampling-creator (sv :sample-type slot-name))
				rng (sv :sample-sigma slot-name) object slot-name))
	(setf (aref sampler-array i) sampler)
	(finally (return sampler-array))))))



(defun --init--sampling-functions (model-object)
  (let+ (((&slots model-parameters-to-marginalize
		  sample-new-parameters) model-object)
	 (no-samplers (length model-parameters-to-marginalize))
	 (array-sampling-fun (--init--sampling-fn-array model-object model-parameters-to-marginalize)))
    (labels ((sample-new-params ()
	       (iter
    		 (for i from 0 below no-samplers)
    		 (funcall (aref array-sampling-fun i))
    		 (finally (return t)))))
      (setf sample-new-parameters #'sample-new-params))))



