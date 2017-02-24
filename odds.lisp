(in-package #:bayesian-analysis)



(defgeneric calculate-odds-ratio-1/2 (model-1 model-2 data &key))


(defmethod calculate-odds-ratio-1/2 ((model-1 model) (model-2 model) (data data) &key)
  (labels ((log-of-prior*likelihood (model)
	     (let+ ((likelihood (initialize-likelihood model data))
		    ((&slots model-prior log-of-all-priors) model)
		    ((&slots constant/log-of-likelihood varying/log-of-likelihood) likelihood)
		    (m-p (log (funcall model-prior)))
		    (l-a-p (funcall log-of-all-priors))
		    (c-l (funcall constant/log-of-likelihood))
		    (v-l (funcall varying/log-of-likelihood)))
	       (+ m-p l-a-p c-l v-l))))
    (let+ ((l1 (log-of-prior*likelihood model-1))
	   (l2 (log-of-prior*likelihood model-2)))
      (exp (- l1 l2)))))


