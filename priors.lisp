(in-package #:bayesian-analysis)



(defun jeffreys-log-lambda (min max object slot-name)
  (let ((ln-max/min (coerce (log (/ max min)) 'double-float)))
    #'(lambda ()
	(log (/ 1 (* (slot-value object slot-name) ln-max/min))))))


(defun jeffreys (min max object slot-name)
  (let ((ln-max/min (coerce (log (/ max min)) 'double-float)))
    (declare (type double-float ln-max/min))
    (/ 1 (* (slot-value object slot-name) ln-max/min))))

(defun uniform-log-lambda (min max object slot-name)
  (declare (ignore object slot-name))
  (let ((1/delta-x (coerce (/ 1 (- max min)) 'double-float)))
    #'(lambda ()
	(declare (type double-float 1/delta-x))
	(log 1/delta-x))))

(defun uniform (min max object slot-name)
  (declare (ignore object slot-name))
  (coerce (/ 1 (- max min)) 'double-float))


(defun make-prior-in-range (min max object slot-name)
  #'(lambda ()
      (declare (type double-float min max))
      (<= min (slot-value object slot-name) max)))



(defparameter *prior-types*
  '((:jeffreys jeffreys-log-lambda :variable jeffreys)
    (:uniform uniform-log-lambda :constant uniform)))


(define-condition unknown-prior-type (error)
  ((type-not-known :accessor type-not-known :initarg :type-not-known
		   :initform :unspecified)))


(defun find-prior-type (type &key (type-list *prior-types*))
  (alexandria:if-let (type (find type type-list :key #'first))
    type (error 'unknown-prior-type :type-not-known type)))

(defun prior-is-constant (type &key (type-list *prior-types*))
  (equal :constant (third (find-prior-type type :type-list type-list))))

(defun get-prior-lambda (type &key (type-list *prior-types*))
  (second (find-prior-type type :type-list type-list)))

(defun get-prior-function (type &key (type-list *prior-types*))
  (fourth (find-prior-type type :type-list type-list)))




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
	(setf (aref prior-in-range-array i)
	      (make-prior-in-range (sv :min slot-name) (sv :max slot-name) object slot-name))
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
		  model-parameters-to-marginalize model-prior) model-object)
	 (no-priors (length model-parameters-to-marginalize))
	 (log-of-constant-priors (--init--get-log-constant-prior model-object model-parameters-to-marginalize))
	 ((&values array-of-log-of-varying-priros no-varying-priors)
	  (--init--prior-array model-object model-parameters-to-marginalize))
	 (array-of-prior-in-range (--init--prior-in-range-array model-object model-parameters-to-marginalize)))
    (labels ((priors-in-range ()
    	       (iter
    		 (for i from 0 below no-priors)
    		 (if (not (funcall (aref array-of-prior-in-range i)))
    		     (return nil))
    		 (finally (return t))))
    	     (log-of-varying-prior ()
    	       (iter
    		 (for i from 0 below no-varying-priors)
    		 (iter:multiply (funcall (aref array-of-log-of-varying-priros i))
				into log-of-p)
    		 (finally (return log-of-p))))
	     (model-prior ()
	       (warn "model-prior Not implemented yet.")
	       1d0))
      (setf priors-in-range #'priors-in-range
	    model-prior #'model-prior
	    constant/log-of-priors (constantly log-of-constant-priors)
	    varying/log-of-priors #'log-of-varying-prior))))
