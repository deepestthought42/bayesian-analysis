(in-package #:bayesian-analysis)


(defgeneric calculate-odds-ratio-2/1 (result-1 result-2 &key prior-result-1 prior-result-2))

(define-condition odds-different-data (error) ()
  (:report (lambda (c stream)
	     (declare (ignore c))
	     (format stream "The two results given had not been
	     initialized with the same data object."))))




(defmethod calculate-odds-ratio-2/1 ((result-1 nlopt-result) (result-2 nlopt-result) 
				     &key (prior-result-1 (constantly 1d0)) (prior-result-2 (constantly 1d0)))
  (if (not (eql (data result-1) (data result-2)))
      (error 'odds-different-data))
  (let+ ((p-m1 (funcall prior-result-1 result-1))
	 (p-m2 (funcall prior-result-2 result-2))
	 (bayes-factor (exp (- (laplacian-approximation-likelihood result-2)
			       (laplacian-approximation-likelihood result-1))))	 )
    (* (/ p-m2 p-m1) bayes-factor)))


