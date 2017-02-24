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
   (log-of-all-priors :accessor log-of-all-priors :initarg :log-of-all-priors
		      :initform (constantly 0d0))
   (sample-new-parameters :initarg :sample-new-parameters :accessor sample-new-parameters)
   (object-documentation :initarg :object-documentation :accessor object-documentation)
   (rng :accessor rng :initarg :rng
	:initform (gsl-cffi:get-random-number-generator gsl-cffi::*mt19937_1999*))
   (model-function :reader model-function)
   (initargs :reader initargs)))



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
    (:sample-type -sample-type)
    (:bin-width -bin-width)
    (:description -description)))



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
						     (bin-width 0.1d0)
						     (description "")
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
      ,(standard-slot :sample-sigma sample-sigma t)
      ,(standard-slot :bin-width bin-width t)
      ,(standard-slot :description description))))



(defun get-all-params (given-param-initializers)
  (let+ ((no-params (length given-param-initializers))
	 (names (mapcar #'first given-param-initializers)))
    (values names no-params)))




(defun make-initialize-after-code (class-name model-function-name parameters documentation model-prior-code)
  `(defmethod initialize-instance :after ((object ,class-name) &rest initargs &key &allow-other-keys)
     (labels ((l-model-prior () ,(if model-prior-code model-prior-code '1d0)))
       (setf (slot-value object 'all-model-parameters) ',parameters
	     (slot-value object 'model-parameters-to-marginalize) (--init--params-to-marginalize object ',parameters)
	     (slot-value object 'model-prior) #'l-model-prior
	     (slot-value object 'initargs) initargs
	     (slot-value object 'model-function) #',model-function-name 
	     (object-documentation object) ,documentation))
     (iter:iter
       (iter:for p in ',parameters)
       (setf (slot-value object p) (coerce (slot-value object p)
					   'double-float)))
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
		 ;; save and restore the old params
		 (setf (aref param-value-array i current-index)
		       (aref param-value-array i (1- current-index))
		       (slot-value model (aref symbol-array i))
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


(defun make-model-class-and-coby-object (name parameters)
  (let+ ((all-slot-specifiers (mapcan #'(lambda (s) (apply #'make-slot-specifiers-for-parameter s))
				      parameters)))
    `((defclass ,name (bayesian-analysis:model)
	(,@all-slot-specifiers))
      (defmethod copy-object ((object ,name))
	(let* ((new-object (apply #'make-instance ',name (initargs object))))
	  (iter:iter
	    (for s in ',(mapcar #'first parameters))
	    (setf (slot-value new-object s)
		  (slot-value object s)))
	  new-object)))))


(defun model-function-name (model-name)
  (alexandria:symbolicate model-name '-model-function))

(defun make-model-function (model-function-name independent-parameters
			    model-parameters body)
  (alexandria:with-gensyms (model-object)
    `(defun ,model-function-name (,@independent-parameters ,model-object)
       (let-plus:let+ (((let-plus:&slots ,@model-parameters) ,model-object))
	 (progn ,@body)))))




(defmacro define-bayesian-model ((name data-type &key model-prior-code documentation) 
				 (&rest model-parameters)
				 (likelihood-type &key constant-likelihood-code varying-likelihood-code)
				 ((&rest independent-parameters) &body body))
  (let+ ((model-function-name (model-function-name name)))
    `(progn
       ,@(make-model-class-and-coby-object name model-parameters)
       ,(make-initialize-after-code name model-function-name
				    (mapcar #'first model-parameters)
				    documentation
				    model-prior-code)
       ,(make-accumulator-method name)
       ,(make-model-function model-function-name independent-parameters
			     (mapcar #'first model-parameters)
			     body)
       ,(make-likelihood-initializer name likelihood-type
				     model-function-name 
				     (bayesian-analysis:get-independent-parameters data-type)
				     (bayesian-analysis:get-dependent-parameters data-type)
				     (bayesian-analysis:get-error-parameters data-type)
				     data-type))))


 

