(in-package #:bayesian-analysis)



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

