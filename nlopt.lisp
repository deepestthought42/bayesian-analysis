(in-package #:bayesian-analysis)

(define-condition nlopt-couldnt-optimize ()
  ((result :accessor nlopt-result :initarg :nlopt-result
		 :initform 'nlopt:+nlopt_failure+)
   (explanation :accessor explanation :initarg :explanation
		:initform (nlopt:docstring 'nlopt:nlopt-result nlopt:+nlopt_failure+)))
  (:report (lambda (c stream)
	    (let+ (((&slots result explanation) c))
	      (format stream "NLopt optimization aborted with errorcode: ~a (~a) "
		      result explanation)))))






(defclass nlopt (algorithm nlopt:config) ())


(defclass nlopt-result (optimization-result)
  ((nlopt-result :initarg :nlopt-result :accessor nlopt-result 
		 :initform (error "Must initialize nlopt-result."))
   (f-val :initarg :f-val :accessor f-val 
	  :initform (error "Must initialize f-val."))
   (likelihood :initarg :likelihood :accessor likelihood 
	       :initform (error "Must initialize likelihood."))
   (model :initarg :model :accessor model 
	  :initform (error "Must initialize model."))))


;; implementation of cffi-nlopt api

(defmethod nlopt:no-dimensions ((likelihood likelihood))
  (length (model-parameters-to-marginalize (model likelihood))))

(defun get-ln-of-p*L (model data)
  (let+ ((likelihood (ba:initialize-likelihood model data))
	 ((&slots varying/log-of-likelihood
		  constant/log-of-likelihood model) likelihood)
	 ((&slots log-of-all-priors model-parameters-to-marginalize) model))
    (lambda ()
      (+
       (funcall varying/log-of-likelihood)
       (funcall constant/log-of-likelihood)
       (funcall log-of-all-priors)))))

(defmethod nlopt:function-to-optimize ((likelihood likelihood))
  (let+ (((&slots model data) likelihood)
	 (ln-of-p*L (get-ln-of-p*L model data))
	 ((&slots model-parameters-to-marginalize) model))
    (labels ((fun (xa)
	       (iter
		 (for p in model-parameters-to-marginalize)
		 (for i initially 0 then (1+ i))
		 (setf (slot-value model p) (cffi:mem-aref xa :double i)))
	       (funcall ln-of-p*L)))
      #'fun)))

(defmethod nlopt:upper-bounds ((likelihood likelihood))
  (let+ (((&slots model) likelihood))
    (iter
      (for p in (model-parameters-to-marginalize model))
      (collect (coerce (slot-value model (w/suffix-slot-category :max p))
		       'double-float)))))

(defmethod nlopt:lower-bounds ((likelihood likelihood))
  (let+ (((&slots model) likelihood))
    (iter
      (for p in (model-parameters-to-marginalize model))
      (collect (coerce (slot-value model (w/suffix-slot-category :min p)) 'double-float)))))

(defmethod nlopt:initial-guess ((likelihood likelihood))
  (let+ (((&slots model) likelihood))
    (iter
      (for p in (model-parameters-to-marginalize model))
      (collect (coerce (slot-value model p) 'double-float)))))


(defmethod find-optimum ((algorithm nlopt) input-model data &key)
  (let+ ((model (copy-object input-model))
	 ((&slots log-of-all-priors) model)
	 (likelihood (ba:initialize-likelihood model data))
	 ((&values retval &ign opt-val) (nlopt:optimization likelihood algorithm)))
    ;; fixme: this needs to go into nlopt
    (if (< retval 0)
	(error 'nlopt-couldnt-optimize
	       :explanation (nlopt:docstring 'nlopt:nlopt-result retval)
	       :nlopt-result (nlopt:get-symbol 'nlopt:nlopt-result retval)))
    (make-instance 'nlopt-result
		   :f-val opt-val
		   :model model
		   :input-model input-model
		   :nlopt-result retval
		   :likelihood likelihood
		   :data data
		   :algorithm algorithm)))


(defun %do-f-max (model-to-optimize data algorithm)
  (let+ ((likelihood (ba:initialize-likelihood model-to-optimize data)))
    (nlopt:optimization likelihood algorithm)))


(defun %laplacian-max-phi (x parameter algorithm model-to-optimize data)
  (setf (slot-value model-to-optimize parameter) x)
  (let+ (((&values &ign &ign f) (%do-f-max model-to-optimize data algorithm))
	 (model-for-hessian (copy-object model-to-optimize))
	 (fun-for-hessian (get-ln-of-p*l model-for-hessian data))
	 (hessian (hessian fun-for-hessian model-for-hessian
			   (get-optimal-delta model-for-hessian) -1d0))
	 (determinant (lla:det hessian))
	 (f (funcall (get-ln-of-p*l model-to-optimize data)))
	 (d (abs (sqrt determinant))))
    (- f (log d))))


(defun laplacian-approximation (nlopt-result parameter no-bins &key on-center)
  (let+ ((marginalize-keyword (alexandria:make-keyword (w/suffix-slot-category :marginalize parameter)))
	 ((&slots model data algorithm) nlopt-result)
	 (model-to-optimize (copy-object model marginalize-keyword nil))
	 (min (coerce (slot-value model (w/suffix-slot-category :min parameter)) 'double-float))
	 (max (coerce (slot-value model (w/suffix-slot-category :max parameter)) 'double-float)))
    (iter
      (with diff = (- max min))
      (for x from min to max by (/ diff no-bins))
      (for ln-of-p*L = (%laplacian-max-phi x parameter algorithm
					   model-to-optimize data))
      (collect x into xs)
      (collect ln-of-p*L into f-of-xs)
      (maximize ln-of-p*L into max-ln-of-p*L)
      (minimize ln-of-p*L into min-ln-of-p*L)
      (finally
       (return
	 (iter
	   (with integral = (sumlogexp f-of-xs))
	   (with result-array = (make-array (length f-of-xs)))
	   (for i initially 0 then (1+ i))
	   (for x in xs)
	   (for f in f-of-xs)
	   (setf (aref result-array i)
		 (list (if on-center (- x (/ (+ max min) 2d0)) x)
		       (exp (- f integral))))
	   (finally (return (values result-array min max)))))))))




(defmethod bin-parameter-values ((result nlopt-result) parameter
				 &key (no-bins 100) (start 0) end (confidence-level 0.69))
  (declare (ignore start end))
  
  (let+ (((&values binned min max)
	  (laplacian-approximation result parameter no-bins)))
    ;; processing
    (iter
      (with median = (/ (- max min) 2))
      (with first-time = t)
      (with median-index = 0)
      (with counts = 0)
      (for (x c) in-sequence binned with-index i)
      (maximize c into max-counts)
      (incf counts c)
      (if (and first-time (>= counts 0.5d0))
	  (setf first-time nil
		median (car (aref binned i)) 
		median-index i))
      (finally
       (let+ (((&values min max)
	       (%calculate-confidance binned 1d0 confidence-level)))
	 (return (values (map 'list #'identity binned)
			 median min max max-counts)))))))

(defmethod get-parameter-results ((result nlopt-result)
				  &key (start 0) end (confidence-level 0.69)
				       (no-bins 50))
  (let+ (((&slots nlopt-result input-model model data) result)
	 ((&slots model-parameters-to-marginalize) input-model)
	 (param-infos (iter
			(for p in model-parameters-to-marginalize)
			(let+ (((&values binned-data median min max max-counts)
				(bin-parameter-values result p
						      :no-bins no-bins
						      :start start :end end
						      :confidence-level confidence-level)))
			  (setf (slot-value model p) (coerce median 'double-float))
			  (collect (make-instance 'parameter-distribution
						  :name p
						  :median median
						  :confidence-level confidence-level
						  :confidence-min min
						  :confidence-max max
						  :max-counts max-counts
						  :binned-data binned-data))))))
    (make-instance 'optimized-parameters
		   :algorithm-result result
		   :parameter-infos param-infos
		   :data data
		   :model model)))
