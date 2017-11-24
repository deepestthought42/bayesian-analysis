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


(defmethod calculate-odds-ratio-1/2 ((model-1 model) (model-2 model) (data data)
				     &key (model-1-prior (constantly 1d0))
					  (model-2-prior (constantly 1d0)))
  (labels ((log-of-prior*likelihood (model model-prior)
	     (let+ ((likelihood (initialize-likelihood model data))
		    ((&slots log-of-all-priors) model)
		    ((&slots constant/log-of-likelihood varying/log-of-likelihood) likelihood)
		    (m-p (log (funcall model-prior)))
		    (l-a-p (funcall log-of-all-priors))
		    (c-l (funcall constant/log-of-likelihood))
		    (v-l (funcall varying/log-of-likelihood)))
	       (+ m-p l-a-p c-l v-l))))
    (let+ ((l1 (log-of-prior*likelihood model-1 model-1-prior))
	   (l2 (log-of-prior*likelihood model-2 model-2-prior)))
      (exp (+ (- l1 l2))))))



