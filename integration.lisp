(in-package #:bayesian-analysis)



(defclass integration () ())



(defgeneric marginalize-parameter (integration model data &key))




(defun get-integration-fun-for-parameter (parameter start end model close-over integration-function-creator)
  ;;; this needs to be documented to be understandable next time someone looks at it
  ;; so, this function will return a function 
  (let+ ((prior (get-prior-for-parameter model parameter))
	 (integrator (funcall integration-function-creator)))
    (labels ((integrand (x)
	       (setf (slot-value model parameter) x)
	       (* (exp (funcall prior))
		  (funcall close-over)))
	     (fun ()
	       (let ((original (slot-value model parameter))
		     (result (funcall integrator #'integrand start end)))
		 (setf (slot-value model parameter) original)
		 result)))
      #'fun)))


(defun integrate-over (model data parameters
		       &key (integration-function-creator #'gsl-cffi:create-integration-function))
  (let+ ((likelihood (ba:initialize-likelihood model data))
	 ((&slots constant/log-of-likelihood varying/log-of-likelihood) likelihood))
    (iter
      (for (p start end) in parameters)
      (for fun initially #'(lambda () (exp (funcall varying/log-of-likelihood)))
	   then (get-integration-fun-for-parameter p start end model fun integration-function-creator))
      (finally
       (return (* (exp (funcall constant/log-of-likelihood))
		  (funcall fun)))))))



