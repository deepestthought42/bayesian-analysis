(in-package #:bayesian-analysis)



(defclass integration () ())



(defgeneric marginalize-parameter (integration model data &key))




(defun get-integration-fun-for-parameter (parameter start end model close-over integration-function-creator)
  ;;; this needs to be documented to be understandable next time someone looks at it
  ;; so, this function will return a function 
  (let+ ((integrator (funcall integration-function-creator)))
    (labels ((integrand (x)
	       (setf (slot-value model parameter) x)
	       (funcall close-over))
	     (fun ()
	       (funcall integrator #'integrand start end)))
      #'fun)))




(defun get-ln-of-joint-prior-lambda (parameters model)
  (let+ ((no-params (length parameters))
	 (arr (make-array no-params :element-type '(function () double-float)
				    :initial-element (constantly 0d0))))      
    (iter
      (for (p) in-sequence parameters with-index i)
      (setf (aref arr i)
	    (get-prior-for-parameter model p)))
    #'(lambda ()
	(iter
	  (for i from 0 below no-params)
	  (sum (funcall (aref arr i)))))))




(defun integrate-over (model data parameters
		       &key (integration-function-creator #'gsl-cffi:create-integration-function))
  (let+ ((model-copy (apply #'ba:copy-object model
			    (iter
			      (for (p start end) in parameters)
			      (collect (alexandria:make-keyword
					(get-slot-name/cat :min p)))
			      (collect start)
			      (collect (alexandria:make-keyword
					(get-slot-name/cat :max p)))
			      (collect end))))
	 (likelihood (ba:initialize-likelihood model-copy data))
	 (ln-of-joint-prior (get-ln-of-joint-prior-lambda parameters model-copy))
	 ((&slots constant/log-of-likelihood varying/log-of-likelihood) likelihood))
    (iter
      (for (p start end) in parameters)
      (for fun initially #'(lambda () (exp
				  (+ (funcall ln-of-joint-prior)
				     (funcall varying/log-of-likelihood))))
	   then (get-integration-fun-for-parameter p start end model-copy fun integration-function-creator))
      (finally
       (return (* (exp (funcall constant/log-of-likelihood))
		  (funcall fun)))))))





