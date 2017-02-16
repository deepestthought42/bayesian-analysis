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

