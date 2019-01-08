(in-package #:bayesian-analysis)

(define-condition odds-different-data (error) ()
  (:report (lambda (c stream)
	     (declare (ignore c))
	     (format stream "The two results given had not been
	     initialized with the same data object."))))


(defgeneric odds-ratio-1/2 (a params-a b params-b))
(defgeneric calculate-marginal-posterior (for-object data params-and-range))


(defmethod calculate-marginal-posterior ((model model) (data data) params-and-range)
  "Given a parameter-result in PARAMETER-RESULT, a set of data that has "
  (let+ ((params (mapcar #'(lambda (p)
			     (let+ (((parameter-name range-start &optional range-end) p)
				    (val (slot-value model parameter-name)))
			       (if range-end 
				   (list parameter-name range-start range-end)
				   (list parameter-name (- val range-start) (+ val range-start)))))
			 params-and-range)))
    (ba::integrate-over model data params)))

(defmethod calculate-marginal-posterior ((parameter-result optimized-parameters) (data data) params-and-range)
  "Given a parameter-result in PARAMETER-RESULT, a set of data that has "
  (let+ ((model (ba:model parameter-result)))
    (calculate-marginal-posterior model data params-and-range)))


(defmethod odds-ratio-1/2 ((result-1 optimized-parameters) params-1
			   (result-2 optimized-parameters) params-2)
  "Given parameter results in RESULT-1 and RESULT-2, this function
will calculate the odds ratio for the result-models of RESULT-1 and
RESULT-2. The marginal posterior will be calculated by numerically
integrating over the params given in PARAMS-1 and PARAMS-2
respectively. PARAMS-1, PARAMS-2 are list of the following form:

{PARAMETER-NAME, DEVIATION}*

where PARAMETER-NAME is the name of the slot in the model to
marginalize over and DEVIATION is half the range, centered on the
optimal value, that will be used to integrate over."
  (let+ ((data-1 (data result-1))
	 (data-2 (data result-2)))
    ;; some sanity checks first
    (if (not (eql data-1 data-2))
	(error 'odds-different-data))
    (/ (calculate-marginal-posterior result-1 data-1 params-1)
       (calculate-marginal-posterior result-2 data-1 params-2))))





