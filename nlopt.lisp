(in-package #:bayesian-analysis)



(defclass nlopt (algorithm nlopt:config) ())


(defclass nlopt-result (optimization-result)
  ((model :initarg :model :accessor model 
	  :initform (error "Must initialize result-model."))
   (nlopt-result :initarg :nlopt-result :accessor nlopt-result 
		 :initform (error "Must initialize nlopt-result."))))


;; implementation of cffi-nlopt api

(defmethod nlopt:no-dimensions ((likelihood likelihood))
  (length (model-parameters-to-marginalize (model likelihood))))

(defmethod nlopt:function-to-optimize ((likelihood likelihood))
  (let+ (((&slots varying/log-of-likelihood
		  constant/log-of-likelihood
		  model) likelihood)
	 ((&slots log-of-all-priors model-parameters-to-marginalize) model))
    (labels ((fun (xa)
	       (iter
		 (for p in model-parameters-to-marginalize)
		 (for i initially 0 then (1+ i))
		 (setf (slot-value model p) (cffi:mem-aref xa :double i)))
	       (+
		(funcall varying/log-of-likelihood)
		(funcall constant/log-of-likelihood)
		(funcall log-of-all-priors))))
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
	 (retval (nlopt:optimization likelihood algorithm)))
    (make-instance 'nlopt-result
		   :model model
		   :input-model input-model
		   :nlopt-result retval
		   :data data
		   :algorithm algorithm)))


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
