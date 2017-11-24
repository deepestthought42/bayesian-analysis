(in-package #:bayesian-analysis)


(declaim (optimize (debug 3) (space 0) (safety 1) (speed 3)))



(defclass metropolis-hastings (algorithm)
  ((no-iterations :accessor no-iterations
		  :initarg :no-iterations :initform 1000)))


(defclass mcmc-optimization-result (optimization-result)
  ((no-accepted-iterations :initarg :no-accepted-iterations :accessor no-accepted-iterations 
			   :initform 0)
   (no-priors-not-in-range :accessor no-priors-not-in-range :initarg :no-priors-not-in-range :initform 0)
   (no-iterations :initarg :no-iterations :accessor no-iterations :initform 0)
   (iteration-accumulator :initarg :iteration-accumulator :accessor iteration-accumulator 
			  :initform (error "Must initialize iteration-accumulator."))))

(defclass mcmc-accumulator ()
  ((parameter-array :initarg :parameter-array :accessor parameter-array 
		    :initform (error "Must initialize parameter-array."))
   (no-iterations :initarg :no-iterations :accessor no-iterations 
		  :initform (error "Must initialize no-iterations."))
   (marginalized-parameters :initarg :marginalized-parameters :accessor marginalized-parameters 
			    :initform (error "Must initialize marginalized-parameters."))
   (no-marginalized-parameters :accessor no-marginalized-parameters :initarg :no-marginalized-parameters)
   (save-current-parameters :initarg :save-current-parameters :accessor save-current-parameters 
			    :initform (error "Must initialize save-current-parameters."))
   (save-last-parameters :initarg :save-last-parameters :accessor save-last-parameters 
			 :initform (error "Must initialize save-last-parameters."))))

(defmethod initialize-instance :after ((object mcmc-accumulator) &key)
  (let+ (((&slots no-marginalized-parameters marginalized-parameters) object))
    (setf no-marginalized-parameters (length marginalized-parameters))))

(defmethod bin-parameter-values ((result mcmc-optimization-result) parameter
				 &key (no-bins 100) (start 0) end (confidence-level 0.69)
				      (fn-on-x #'identity))
  (let+ (((&slots iteration-accumulator) result)
	 ((&slots marginalized-parameters parameter-array
		  no-iterations) iteration-accumulator)
	 (end (if (not end) no-iterations end))
	 (no-data-points (- no-iterations start))
	 (pos (position parameter marginalized-parameters))
	 ((&values bin-width no-bins) (%get-bin-width/no-bins parameter-array pos start end no-bins)))
    (if (not pos)
	(error "Parameter: ~a was not marginalized over." parameter))
    (let+ (((&values vals min max)
	    (iter
	      (for i from start below end)
	      (for val = (* bin-width (floor (aref parameter-array pos i) bin-width)))
	      (maximize val into max)
	      (minimize val into min)
	      (collect val into vals)
	      (finally (return (values vals min max)))))
	   (binned (make-array (1+ no-bins)
				  :initial-contents
				  (iter
				    (for i from 0 to no-bins)
				    (for x = (+ min (* (+ i 0.5) bin-width)))
				    (collect (list x 0))))))
      ;; binning
      (iter
	(for v in vals)
	(incf (cadr (aref binned (round (- v min) bin-width)))))
      ;; processing
      (iter
	(with median = (/ (- max min) 2))
	(with first-time = t)
	(with median-index = 0)
	(with counts = 0)
	(for (x c) in-sequence binned with-index i)
	(maximize c into max-counts)
	(incf counts c)
	(if (and first-time (>= counts (/ no-data-points 2)))
	    (setf first-time nil
		  median (car (aref binned i))
		  median-index i))
	(finally
	 (let+ (((&values min max)
		 (%calculate-confidance binned no-data-points confidence-level)))
	   (return (values (map 'list #'(lambda (n) (list (funcall fn-on-x (first n))
						     (/ (second n) no-data-points)))
				binned)
			   median min max (/ max-counts no-data-points)))))))))


(defmethod get-parameter-results ((result mcmc-optimization-result)
				  &key (start 0) end (confidence-level 0.9545)
				       (no-bins 50) (fn-on-x #'identity))
  (let+ (((&slots iteration-accumulator input-model no-iterations data) result)
	 ((&slots model-parameters-to-marginalize) input-model)
	 (model (copy-object input-model))
	 (end (if end end no-iterations))
	 (param-infos (iter
			(for p in model-parameters-to-marginalize)
			(for param-dist = (make-parameter-distribution result p no-bins start end
								       confidence-level no-iterations
								       :fn-on-x fn-on-x))
			(setf (slot-value model p) (coerce (median param-dist) 'double-float))
			(collect param-dist))))
    (make-instance 'optimized-parameters
		   :algorithm-result result
		   :parameter-infos param-infos
		   :data data
		   :model model)))



(defmethod find-optimum ((algorithm metropolis-hastings) input-model data
			 &key random-numbers (force-accept-after 100))
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
	(declare (type fixnum accepted-iterations i no-iterations discarded-priors) 
		 (type (simple-array (double-float 0d0)) random-numbers)
		 (type double-float log-mh-ratio log-u
		       nominator denominator)
		 (ftype (function ((function () double-float)) double-float) ->)
		 (ftype (function ((function () t)) t) !>))
	(with accepted-iterations = 0)
	(with discarded-priors = 0)
	(with denominator = (+ (-> varying/log-of-priors)
			       (-> varying/log-of-likelihood)))
	(for i from 0 below no-iterations)
	(!> sample-new-parameters)
	(when (not (!> priors-in-range))
	  (progn
	    (incf discarded-priors)
	    (!> save-last-parameters))
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
		   (make-instance 'mcmc-optimization-result
				  :algorithm algorithm
				  :no-accepted-iterations accepted-iterations
				  :no-priors-not-in-range discarded-priors
				  :no-iterations no-iterations
				  :data data
				  :input-model input-model
				  :iteration-accumulator accumulator)))))))
 
