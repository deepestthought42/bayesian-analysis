(in-package #:bayesian-analysis)




(defun calculate-marginal-posterior (parameter-result data params-and-range) 
  (let+ ((model (ba:model parameter-result))
	 (params (mapcar #'(lambda (p)
			     (let+ (((parameter-name range/2) p)
				    (val (slot-value model parameter-name)))
			       (list parameter-name (- val range/2) (+ val range/2))))
			 params-and-range)))
    (ba::integrate-over model data params)))

(defun odds-ratio-1/2 (model-1 params-1 model-2 params-2
		       data &key (no-bins 50) (confidence-level 0.9545))
  (labels ((calc (model params)
	     (calculate-marginal-posterior
	      (ba:get-parameter-results model :no-bins no-bins
					      :confidence-level confidence-level)
	      data params)))
    (/ (calc model-1 params-1)
       (calc model-2 params-2))))





