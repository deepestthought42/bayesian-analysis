(in-package #:bayesian-analysis)


;(declaim (optimize (debug 3) (space 0) (safety 1) (speed 3)))

(defclass metropolis-hastings (mcmc-algorithm) ())




(defmethod solve-for-parameters ((algorithm metropolis-hastings) input-model data
				 &key random-numbers )
  (labels ((!> (fun)
	     (declare (type (function () t) fun))
	     (funcall fun))
	   (-> (fun)
	     (declare (type (function () double-float) fun))
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
	(declare (type fixnum accepted-iterations i no-iterations)
		 (type (simple-array double-float) random-numbers)
		 (type double-float log-mh-ratio log-u
		       nominator denominator))
	(with accepted-iterations = 0)
	(with denominator = (+ (-> varying/log-of-priors)
			       (-> varying/log-of-likelihood)))
	(for i from 0 below no-iterations)
	(!> sample-new-parameters)
	(when (not (!> priors-in-range))
	  (!> save-last-parameters)
	  (next-iteration))
	(for nominator = (+ (-> varying/log-of-priors)
			    (-> varying/log-of-likelihood)))
	(for log-mh-ratio = (- nominator denominator))
	(for log-u = (the double-float (log (aref random-numbers i))))
	(if (> log-u log-mh-ratio)
	    (!> save-last-parameters)
	    (progn
	      (incf accepted-iterations)
	      (!> save-current-parameters)
	      (setf denominator nominator)))

	(finally (return
		   (make-instance 'mcmc-parameter-result
				  :algorithm algorithm
				  :no-accepted-iterations accepted-iterations
				  :no-iterations no-iterations
				  :data data
				  :input-model input-model
				  :iteration-accumulator accumulator)))))))
 
