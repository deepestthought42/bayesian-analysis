(in-package #:bayesian-analysis)


(declaim (optimize (debug 3) (speed 1) (space 0)))


(defclass levenberg-marquardt (algorithm)
  ((error-absolute :accessor error-abso :initarg :error-abso :initform 1d-6)
   (error-relative :accessor error-relative :initarg :error-relative :initform 1d-6)
   (max-iterations :accessor max-iterations :initarg :max-iterations :initform 100)))


(defclass levenberg-marquardt-parameter-result (parameter-result)
  ((result-model :initarg :result-model :accessor result-model 
		 :initform (error "Must initialize result-model."))
   (iterations :initarg :iterations :accessor iterations 
	       :initform (error "Must initialize iterations."))
   (gsl-test-delta-status :initarg :gsl-test-delta-status :accessor gsl-test-delta-status 
			  :initform (error "Must initialize gsl-test-delta-status."))
   (gsl-it-status :initarg :gsl-it-status :accessor gsl-it-status 
		  :initform (error "Must initialize gsl-it-status."))))


(defclass iteration ()
  ((no-iteration :accessor no-iteration :initarg :no-iteration :initform 0)
   (chi^2 :initarg :chi^2 :accessor chi^2 
	  :initform (error "Must initialize chi^2."))
   (d-o-f :initarg :d-o-f :accessor d-o-f 
	  :initform (error "Must initialize d-o-f."))
   (model :initarg :model :accessor model 
	  :initform (error "Must initialize model."))
   (covariance-matrix :initarg :covariance-matrix :accessor covariance-matrix 
		      :initform (error "Must initialize covariance-matrix."))))



(defun %analyze-iteration (i no-params solver data-count gsl-vector->params params)
  (let ((covar (gsl-cffi:gsl-matrix-alloc no-params no-params)))
    (unwind-protect
	 (progn
	   (cffi:with-foreign-slots ((gsl-cffi:J gsl-cffi:f gsl-cffi:x)
				     solver
				     (:struct gsl-cffi:fdf-solver))
	     solver
	     (gsl-cffi:gsl-multifit-covar gsl-cffi:J 0d0 covar)
	     (let* ((cov-matrix (make-array (list no-params no-params)
					    :element-type 'double-float))
		    ;; gsl-blas-dnrm2 calculates: \sqrt\sum_i{f_i^2}
		    (chi (gsl-cffi:gsl-blas-dnrm2 gsl-cffi:f))
		    (dof (- data-count no-params)))
	       (iter
		 (for row from 0 below no-params)
		 (iter
		   (for col from 0 below no-params)
		   (setf (aref cov-matrix row col)
			 (gsl-cffi:gsl-matrix-get covar row col))))
	       (make-instance 'iteration
			      :no-iteration i
			      :chi^2 (* chi chi)
			      :d-o-f dof
			      :model (copy-object (funcall gsl-vector->params gsl-cffi:x params))
			      :covariance-matrix cov-matrix))))
      (gsl-cffi:gsl-matrix-free covar))))

(defun %fit-iteration (no-params data-count gsl-vector->params params /solver/
		       error-abs error-rel max-iterations)
  (labels ((analyze-it (i)
	     (%analyze-iteration i no-params
				 /solver/ data-count
				 gsl-vector->params
				 params))
	   (test-delta ()
	     (cffi:with-foreign-slots ((gsl-cffi:dx gsl-cffi:x)
				       /solver/
				       (:struct gsl-cffi:fdf-solver)) 
	       (gsl-cffi:gsl-multifit-test-delta gsl-cffi:dx gsl-cffi:x error-abs error-rel))))
    (iter
      (for i initially 1 then (1+ i))
      (for status-it = (gsl-cffi:gsl-multifit-fdfsolver-iterate /solver/))

      (unless (= gsl-cffi:+GSL_SUCCESS+ status-it)
	(collect (analyze-it i) into it-res)
	(return (values it-res 0 status-it))
	(finish))
    
      (for status = (test-delta))

      (collect (analyze-it i) into it-res)
    
      (while (and (<= i max-iterations)
		  (eql status gsl-cffi:+GSL_CONTINUE+)))
    
      (finally (return (values it-res (if status status 0) status-it))))))

(defun %create-gsl-fit-lambda (gsl-vector->params model data)
  (let+ (((&slots y_i-f_i) model)
	 ((&slots no-data-points) data))
    #'(lambda (x ldata f)
	(declare (ignore ldata))
	(funcall gsl-vector->params x model)
	(iter:iter
	  (iter:for i from 0 below no-data-points)
	  (for Q = (funcall y_i-f_i i model data))
	  (gsl-cffi:gsl-vector-set f i Q))
	;; well, if something goes wrong here we catch it on the lisp
	;; side, so:
	gsl-cffi:+GSL_SUCCESS+)))

(defun make--gsl-vector->params (model)
  (let+ ((slots-to-fit (model-parameters-to-marginalize model)))
    #'(lambda (vec model)
	(iter:iter
	  (iter:for s in slots-to-fit)
	  (iter:for i initially 0 then (1+ i))
	  (setf (slot-value model s)
		(gsl-cffi:gsl-vector-get vec i)))
	model)))

(defun make--initial-params->gsl-vector (model)
  (let ((slots-to-fit (model-parameters-to-marginalize model)))
    #'(lambda (vec model)
	(iter:iter
	  (iter:for s in slots-to-fit)
	  (iter:for i initially 0 then (1+ i))
	  (gsl-cffi:gsl-vector-set vec i (slot-value model s)))
	vec)))




(defun %check-input (model)
  (cond
    ((not (model-parameters-to-marginalize model))
     (error 'no-parameters-to-marginalize))))


(defparameter *cb-fit-lambda*
  #'(lambda (x ldata f)
      (declare (ignore x ldata f))
      (error "Did you forget to register a callback for the fit function ?")
      gsl-cffi:+GSL_FAILURE+))

(cffi:defcallback cb-fit-function
    :int
    ((x gsl-cffi:gsl-vector)
     (ldata :pointer)
     (f gsl-cffi:gsl-vector))
  ;; call thread local variable
  (funcall *cb-fit-lambda* x ldata f))




(defun %levenberg-marquardt-max-likelihood (algorithm model data)
  (gsl-cffi:set-error-handler (cffi:callback gsl-cffi:gsl-error-handler))
  (let+ (((&slots error-absolute error-relative max-iterations) algorithm)
	 (no-fit-params (length (model-parameters-to-marginalize model))) 
	 (gsl-vector->params (make--gsl-vector->params model))
	 (init-params->gsl-vector (make--initial-params->gsl-vector model))
	 ((&slots no-data-points) data)
	 (*cb-fit-lambda* (%create-gsl-fit-lambda gsl-vector->params model data)))
    (gsl-cffi:with-fdf-solver (/solver/ no-data-points no-fit-params
			       #'(lambda (gsl-vector) (funcall init-params->gsl-vector gsl-vector model))
			       (cffi:callback cb-fit-function))
      (multiple-value-bind (fit-iterations gsl-test-delta-status gsl-it-status)
	  (%fit-iteration no-fit-params no-data-points gsl-vector->params model
			  /solver/ error-absolute error-relative max-iterations)
	(let* ((fit-iterations (sort fit-iterations #'> :key #'no-iteration))
	       (res (first fit-iterations))
	       (m (model res)))
	  (make-instance 'levenberg-marquardt-parameter-result
			 :algorithm algorithm
			 :input-model model
			 :data data
			 :result-model m
			 :iterations fit-iterations
			 :gsl-test-delta-status gsl-test-delta-status
			 :gsl-it-status gsl-it-status))))))




(defmethod solve-for-parameters ((algorithm levenberg-marquardt) input-model data &key)
  (%check-input input-model)
  (%levenberg-marquardt-max-likelihood algorithm input-model data))


(defun %laplace-approximate-posterior (model parameter no-bins)
  ""
  (labels ((val (category)
	       (slot-value model (w/suffix-slot-category category parameter))))
    (let+ ((min (val :min))
	   (max (val :max))))))

(defun %make-parameter-result ()
  (let+ (((&values binned-posterior median min max max-counts)))
    
    (make-instance 'solved-parameter-result
		   :name p
		   :median median
		   :confidence-level confidence-level
		   :confidence-min min
		   :confidence-max max
		   :max-counts max-counts
		   :binned-data binned-data)))

(defmethod get-parameter-results ((result levenberg-marquardt-parameter-result)
				  &key (confidence-level 0.69)
				       (no-bins 50))
  (let+ (((&slots input-model result-model data) result)
	 ((&slots model-parameters-to-marginalize) input-model)
	 (param-infos (iter
			(for p in model-parameters-to-marginalize)
			(%make-parameter-result )
			)))
    (make-instance 'solved-parameters
		   :algorithm-result result
		   :parameter-infos param-infos
		   :data data
		   :model model)))
