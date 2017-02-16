(in-package #:bayesian-analysis)


(define-condition initial-parameter-out-of-range ()
  ((parameter :accessor parameter :initarg :parameter :initform :unknown)))


(defclass model ()
  ((all-model-parameters :reader all-model-parameters)
   (model-parameters-to-marginalize :reader model-parameters-to-marginalize)
   (model-prior :accessor model-prior :initarg :model-prior
		:initform (constantly 1d0))
   (constant/log-of-priors :accessor constant/log-of-priors :initarg :constant/log-of-priors
			   :initform (constantly 0d0))
   (varying/log-of-priors :accessor varying/log-of-priors :initarg :varying/log-of-priors
			  :initform (constantly 0d0))
   (priors-in-range :accessor priors-in-range :initarg :priors-in-range
		    :initform (constantly nil))
   (sample-new-parameters :initarg :sample-new-parameters :accessor sample-new-parameters)
   (object-documentation :initarg :object-documentation :accessor object-documentation)
   (rng :accessor rng :initarg :rng
	:initform (gsl-cffi:get-random-number-generator gsl-cffi::*mt19937_1999*))))



(defgeneric initialize-likelihood (model data))


(defgeneric copy-object (object))
(defgeneric cleanup-object (object))
(defgeneric get-slot-value (model slot))





(defun %s (things &rest other-things)
  (apply #'alexandria:symbolicate things other-things))

(defparameter *slot-type<->suffixes*
  '((:marginalize -marginalize)
    (:prior-type -prior-type)
    (:min -min)
    (:max -max)
    (:sample-sigma -sample-sigma)
    (:sample-type -sample-type)))



(defun w/suffix-slot-category (category name &key is-key
					(list-of-suffixes *slot-type<->suffixes*))
  (alexandria:if-let (suffix (find category list-of-suffixes :key #'first))
    (let ((name-with-suffix (alexandria:symbolicate name (second suffix))))
      (if is-key
	  (alexandria:make-keyword name-with-suffix)
	  name-with-suffix))))



(defun make-slot-specifiers-for-parameter (name &key default min max
						     (sample-sigma 1d0)
						     (sample-type :gaussian)
						     (marginalize t)
						     (prior-type :jeffreys))
  (labels ((make-default-param (name default-value &optional coerce)
	     (if (not default-value)
		 `(error ,(format nil "No value known for parameter: ~a" name))
		 (if coerce
		     `,(coerce default-value 'double-float)
		     `,default-value)))
	   (n (cat &key is-key)
	     (w/suffix-slot-category cat name :is-key is-key))
	   (standard-slot (category default &optional coerce)
	     (let ((slot-name (n category)))
	       `(,slot-name
		 :reader slot-name
		 :initarg ,(n category :is-key t)
		 :initform ,(make-default-param slot-name default coerce)))))
    `((,name :accessor ,name :initarg ,(alexandria:make-keyword name)
	     :initform ,(make-default-param name default t))
      ,(standard-slot :marginalize marginalize)
      ,(standard-slot :prior-type prior-type)
      ,(standard-slot :min min t)
      ,(standard-slot :max max t)
      ,(standard-slot :sample-type sample-type)
      ,(standard-slot :sample-sigma sample-sigma t))))



(defun get-all-params (given-param-initializers)
  (let+ ((no-params (length given-param-initializers))
	 (names (mapcar #'first given-param-initializers)))
    (values names no-params)))




(defun make-initialize-after-code (class-name parameters documentation)
  `(defmethod initialize-instance :after ((object ,class-name) &key)
     (setf (slot-value object 'all-model-parameters) ',parameters
	   (slot-value object 'model-parameters-to-marginalize)
	   (--init--params-to-marginalize object ',parameters)
	   (object-documentation object) ,documentation)
     (--init--priors object)
     (--init--sampling-functions object)))


(defun --acc--initialize (model no-iterations)
  (let+ ((no-iterations (1+ no-iterations))
	 ((&slots model-parameters-to-marginalize) model)
	 (no-params (length model-parameters-to-marginalize))
	 (param-value-array  (make-array (list no-params no-iterations)
				  :initial-element 0d0
				  :element-type 'double-float))
	 (symbol-array (make-array no-params :element-type 'symbol
				   :initial-contents model-parameters-to-marginalize))
	 (current-index 0))
    (labels ((save-current ()
	       (iter
		 (for i from 0 below no-params)
		 (setf (aref param-value-array i current-index)
		       (slot-value model (aref symbol-array i)))
		 (finally
		  (incf current-index)
		  (return param-value-array))))
	     (save-last ()
	       (iter
		 (for i from 0 below no-params)
		 (setf (aref param-value-array i current-index)
		       (aref param-value-array i (1- current-index)))
		 (finally
		  (incf current-index)
		  (return param-value-array)))))
      ;; save initial state and get limits strictly correct
      (save-current)
      (make-instance 'mcmc-accumulator
		     :save-last-parameters #'save-last
		     :save-current-parameters #'save-current
		     :no-iterations no-iterations
		     :no-marginalized-parameters no-params
		     :marginalized-parameters model-parameters-to-marginalize
		     :parameter-array param-value-array))))


(defun make-accumulator-method (name)
  `(defmethod bayesian-analysis:initialize-accumulator ((object ,name) no-iterations)
     (--acc--initialize object no-iterations)))


(defun make-model-class (name parameters)
  `(defclass ,name (bayesian-analysis:model)
     (,@(mapcan #'(lambda (s) (apply #'make-slot-specifiers-for-parameter s)) parameters))))



(defun make-likelihood-initializer (model fn-param-list function-body)
  `(defmethod bayesian-analysis:initialize-likelihood ((model ,name) (data bayesian-analysis:data))))

(defmacro define-bayesian-model ((name &optional documentation) 
				 (&rest model-parameters)
				 (log-likelihood)
				 (x-data &rest more-data)
				 &body body)
  `(progn
     ,(make-model-class name model-parameters)
     ,(make-initialize-after-code name (mapcar #'first model-parameters) documentation)
     ,(make-accumulator-method name)))




#+nil
(define-bayesian-model (bayesian-konig)
    ((a :default 2 :min 1 :max 20)
     (b :prior-type :uniform :default 2 :min 1 :max 20))
    ((x))
  )

#+nil
(defparameter *test-model* (make-instance 'bayesian-konig))

(ba:initialize-accumulator *test-model* 1000)

;(funcall (sample-new-parameters *test-model*))



