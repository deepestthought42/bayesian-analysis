(in-package #:bayesian-analysis)



(defclass integration () ())



(defgeneric marginalize-parameter (integration model data &key))




(defun get-integration-fun-for-parameter (parameter model close-over integration-function)
  ;;; this needs to be documented to be understandable next time someone looks at it
  ;; so, this function will return a function 
  (labels ((v (cat)
	     (slot-value model (w/suffix-slot-category cat parameter))))
    (let+ ((start (v :min))
	   (end (v :max))
	   (prior (get-prior-for-parameter model parameter)))
      (labels ((integrand (x)
		 (setf (slot-value model parameter) x)
		 (* (exp (funcall prior))
		    (funcall close-over)))
	       (fun ()
		 (let ((original (slot-value model parameter))
		       (result (funcall integration-function #'integrand start end)))
		   (setf (slot-value model parameter) original)
		   result)))
	#'fun))))

  
(defun integrate-over (model data parameters
		       &key (integration-function (gsl-cffi:create-integration-function)))
  (let+ ((likelihood (ba:initialize-likelihood model data))
	 ((&slots constant/log-of-likelihood varying/log-of-likelihood) likelihood))
    (iter
      (for p in parameters)
      (for fun initially #'(lambda () (* 1d200 (exp (funcall varying/log-of-likelihood))))
	   then (get-integration-fun-for-parameter p model fun integration-function))
      (finally
       (return (* (exp (funcall constant/log-of-likelihood))
		  (funcall fun)))))))



