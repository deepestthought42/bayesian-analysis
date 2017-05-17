(in-package #:bayesian-analysis)



(defclass nlopt (algorithm nlopt:config) ())


(defclass nlopt-result (optimization-result)
  ((model :initarg :model :accessor model 
	  :initform (error "Must initialize result-model."))
   (nlopt-result :initarg :nlopt-result :accessor nlopt-result 
		 :initform (error "Must initialize nlopt-result."))
   (f-val :initarg :f-val :accessor f-val 
	  :initform (error "Must initialize f-val."))
   (likelihood :initarg :likelihood :accessor likelihood 
	       :initform (error "Must initialize likelihood."))))


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


(defun %laplacian-max-phi (x parameter algorithm model-to-optimize data logarithmic)
  (setf (slot-value model-to-optimize parameter) x)
  (let+ (((&values &ign &ign f) (%do-f-max model-to-optimize data algorithm))
	 (model-for-hessian (copy-object model-to-optimize))
	 (fun-for-hessian (get-ln-of-p*l model-for-hessian data))
	 (hessian (hessian fun-for-hessian model-for-hessian
			   (get-optimal-delta model-for-hessian) -1d0))
	 ((&values &ign determinant) (math-utils:invert-matrix hessian))
	 (f (funcall (get-ln-of-p*l model-to-optimize data)))
	 (d (sqrt determinant)))
    (if logarithmic
	(values (- f (log d)) f (log d))
	(values (exp (- f (log d))) (exp f) d))))


(defun laplacian-approximation (nlopt-result parameter no-bins &key (collect-log nil))
  (let+ ((marginalize-keyword (alexandria:make-keyword (w/suffix-slot-category :marginalize parameter)))
	 ((&slots model data algorithm) nlopt-result)
	 (model-to-optimize (copy-object model marginalize-keyword nil)))
    (let+ ((min (coerce (slot-value model (w/suffix-slot-category :min parameter)) 'double-float))
	   (max (coerce (slot-value model (w/suffix-slot-category :max parameter)) 'double-float)))
      (iter
	(with diff = (- max min))
	(for x from min to max by (/ diff no-bins))
	(for (values f-d f d) =
	     (%laplacian-max-phi x parameter algorithm model-to-optimize data collect-log))
	(collect (list x f-d f d))))))




(defmethod get-parameter-results ((result nlopt-result)
				  &key (start 0) end (confidence-level 0.69)
				       (no-bins 50))
  (let+ (((&slots nlopt-result input-model model data) result)
	 ((&slots model-parameters-to-marginalize) input-model)
	 (param-infos (iter
			(for p in model-parameters-to-marginalize)
			(collect (make-instance 'parameter-distribution
						:name p
						:median (coerce (slot-value model p) 'double-float)
						:confidence-level confidence-level
						:confidence-min 0d0
						:confidence-max 0d0
						:max-counts 0
						:binned-data nil)))))
    (make-instance 'optimized-parameters
		   :algorithm-result result
		   :parameter-infos param-infos
		   :data data
		   :model model)))
