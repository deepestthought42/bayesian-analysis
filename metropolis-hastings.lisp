(in-package #:bayesian-analysis)




(defclass metropolis-hastings (mcmc-algorithm) ())




(defmethod solve-for-parameters ((algorithm metropolis-hastings) input-model data
				 &key random-numbers )
  (labels ((-> (fun)
	     (funcall fun)))
    (let+ (((&slots no-iterations) algorithm)
	   (model (copy-object input-model))
	   (random-numbers (if (not random-numbers)
			       (gsl-cffi:get-array-random-uniform no-iterations)))
	   ((&slots varying/log-of-priors
		    priors-in-range
		    sample-new-parameters) model)
	   (likelihood (initialize-likelihood model data))
	   ((&slots varying/log-of-likelihood) likelihood)
	   (accumulator (initialize-accumulator model no-iterations))
	   ((&slots save-current-parameters
		    save-last-parameters) accumulator))
      (declare (type (function () boolean)
		     sample-new-parameters
		     save-current-parameters
		     save-last-parameters)
	       (type (function () double-float)
		     varying/log-of-likelihood
		     varying/log-of-priors))
      (iter
	(with accepted-iterations = 0)
	(with nominator = (* (-> varying/log-of-priors)
			     (-> varying/log-of-likelihood)))

	(for i from 0 below no-iterations)
	(-> sample-new-parameters)
	
	(when (not (-> priors-in-range))
	  (-> save-last-parameters)
	  (next-iteration))
	(for denominator = (* (-> varying/log-of-priors)
			      (-> varying/log-of-likelihood)))
	
	(for log-mh-ratio = (- (- nominator denominator)))
	(for log-u = (log (aref random-numbers i)))
	(if (> log-u log-mh-ratio)
	    (-> save-last-parameters)
	    (progn
	      (incf accepted-iterations)
	      (-> save-current-parameters)
	      (setf nominator denominator)))

	(finally (return
		   (make-instance 'mcmc-parameter-result
				  :algorithm algorithm
				  :no-accepted-iterations accepted-iterations
				  :no-iterations no-iterations
				  :data data
				  :input-model input-model
				  :iteration-accumulator accumulator)))))))
 
