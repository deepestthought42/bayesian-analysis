(in-package #:bayesian-analysis)



(defclass likelihood ()
  ((constant/log-of-likelihood :accessor constant/log-of-likelihood :initarg
			       :constant/log-of-likelihood
			       :initform (constantly 0d0))
   (varying/log-of-likelihood :accessor varying/log-of-likelihood :initarg :varying/log-of-likelihood
			      :initform (constantly 0d0))))

(defgeneric likelihood (likelihood)
  (:method ((l likelihood))
    (exp (+ (funcall (constant/log-of-likelihood l))
	    (funcall (varying/log-of-likelihood l))))))


(defparameter *likelihood-types*
  '((:d_i=f_i+gaussian_error_i)))

(defun make-likelihood-codes (likelihood-type 
			      model-function-name model-object-name
			      data-independent-parameters
			      data-dependent-parameters
			      data-error-parameters
			      data-object-name)
  (case likelihood-type
    (:d_i=f_i+gaussian_error_i
     (if (not (= 1 (length data-dependent-parameters) (length data-error-parameters)))
	 (error "For a simple d_i = f_i + e_i, only one error
	 parameter and one dependent parameter is supported."))
     (let+ ((all-data-parameters (append data-independent-parameters
					 data-dependent-parameters
					 data-error-parameters)))
       (values
	(alexandria:with-gensyms (f_i Q fun)
	  `(labels ((,fun (,@all-data-parameters)
		      (let* ((,f_i (,model-function-name ,@data-independent-parameters
							 ,model-object-name))
			     (,Q (- ,@data-dependent-parameters ,f_i)))
			(expt (/ ,Q ,@data-error-parameters) 2))))
	     (let-plus:let+ (((let-plus:&slots ,@all-data-parameters) ,data-object-name)
			     (d_i (map 'simple-vector #',fun ,@all-data-parameters)))
	       (reduce #'+ d_i))))
	`(progn
	   (warn "constant likelihood not implemented correctly.")
	   (constantly 1d0)))))))



(defun make-likelihood-initializer (model-name likelihood-type
				    model-function-name
				    data-independent-parameters
				    data-dependent-parameters
				    data-error-parameters
				    data-type)
  (alexandria:with-gensyms (data-object-name model-object-name)
    (let+ (((&values varying constant)
	    (make-likelihood-codes likelihood-type 
				  model-function-name model-object-name
				  data-independent-parameters data-dependent-parameters
				  data-error-parameters data-object-name)))
      `(defmethod bayesian-analysis:initialize-likelihood ((,model-object-name ,model-name)
							   (,data-object-name ,data-type))
	 (labels ((varying ()
		    ,varying)
		  (constant ()
		    ,constant))
	   (make-instance 'likelihood
			  :varying/log-of-likelihood #'varying
			  :constant/log-of-likelihood #'constant))))))

