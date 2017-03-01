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


(defun make-f_i+gaussian-error_i-likelihood (model-function-name model-object-name
					     data-independent-parameters
					     data-dependent-parameters
					     data-error-parameters
					     data-object-name)
  (let+ ((all-data-parameters (append data-independent-parameters
				      data-dependent-parameters
				      data-error-parameters)))

    (alexandria:with-gensyms (f_i Q fun N)
      (values
       `(labels ((,fun (,@all-data-parameters)
		   (declare (type (double-float) ,@all-data-parameters))
		   (let* ((,f_i (the double-float
				     (,model-function-name ,@data-independent-parameters
							   ,model-object-name)))
			  (,Q (the double-float
				   (- ,@data-dependent-parameters ,f_i))))
		     (expt (/ ,Q ,@data-error-parameters) 2))))
	  (let-plus:let+ (((let-plus:&slots ,@all-data-parameters) ,data-object-name))
	    (declare (type (simple-array double-float) ,@all-data-parameters))
	    (- (iter:iter
		 (declare (type fixnum i))
		 (iter:for i from 0 below (ba:no-data-points ,data-object-name))
		 (iter:sum (funcall #',fun ,@(mapcar #'(lambda (p) `(aref ,p i)) all-data-parameters)))))))
       `(let-plus:let+ (((let-plus:&slots ,@data-error-parameters) ,data-object-name)
			(,N (coerce (length ,@data-error-parameters) 'double-float)))
	  (declare (type (simple-array double-float) ,@data-error-parameters)
		   (type double-float ,N))
	  (log (/ 1 (* (expt (* 2 pi) (/ ,N 2))
		       (iter:iter
			 (declare (type fixnum i))
			 (iter:for i from 0 below (ba:no-data-points ,data-object-name))
			 (iter:multiply (aref ,@data-error-parameters i)))))))))))


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
     (make-f_i+gaussian-error_i-likelihood model-function-name model-object-name
					   data-independent-parameters
					   data-dependent-parameters
					   data-error-parameters
					   data-object-name))))



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

