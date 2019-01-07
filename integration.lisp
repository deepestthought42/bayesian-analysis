(in-package #:bayesian-analysis)



(defclass integration () ())



(defgeneric marginalize-parameter (integration model data &key))




(defun get-integration-fun-for-parameter (parameter start end model close-over integration-function-creator)
  ;;; this needs to be documented to be understandable next time someone looks at it
  ;; so, this function will return a function 
  (let+ ((integrator (funcall integration-function-creator)))
    (labels ((integrand (x)
	       (setf (slot-value model parameter) x)
	       (funcall close-over))
	     (fun ()
	       (funcall integrator #'integrand start end)))
      #'fun)))




(defun get-ln-of-joint-prior-lambda (parameters model)
  (let+ ((no-params (length parameters))
	 (arr (make-array no-params :element-type '(function () double-float)
				    :initial-element (constantly 0d0))))      
    (iter
      (for (p) in-sequence parameters with-index i)
      (setf (aref arr i)
	    (get-prior-for-parameter model p)))
    #'(lambda () (iter (for i from 0 below no-params) (sum (funcall (aref arr i)))))))




(defun integrate-over (model data parameters
		       &key (integration-function-creator #'gsl-cffi:create-integration-function))
  "fixme: document"
  (let+ ((model-copy (apply #'ba:copy-object model
			    (iter
			      (for (p start end) in parameters)
			      (collect (alexandria:make-keyword
					(get-slot-name/cat :min p)))
			      (collect start)
			      (collect (alexandria:make-keyword
					(get-slot-name/cat :max p)))
			      (collect end))))
	 (likelihood (ba:initialize-likelihood model-copy data))
	 (ln-of-joint-prior (get-ln-of-joint-prior-lambda parameters model-copy))
	 ((&slots constant/log-of-likelihood varying/log-of-likelihood) likelihood))
    (labels ((d (n) (coerce n 'double-float)))
      (iter
	(for fun initially #'(lambda () (exp
				    (+ (funcall ln-of-joint-prior)
				       (funcall varying/log-of-likelihood))))
	     ;; get-integration-fun-for-parameter creates a recursive integration function
	     then (get-integration-fun-for-parameter p (d start) (d end) model-copy fun integration-function-creator))
	(for (p start end) in parameters)
	(finally (return (* (exp (funcall constant/log-of-likelihood)) (funcall fun))))))))





(defun parameter-pdf-integrate (parameter no-bins parameters-to-marginalize model data
				&key (integration-function-creator #'gsl-cffi:create-integration-function)
				     normalize)
  (labels ((d (n) (coerce n 'double-float)))
    (let+ (((p-name begin end) parameter)
	   (new-model (copy-object model))
	   (log-prior (get-prior-for-parameter new-model p-name)))
      (iter
	(for x from (d begin) to (d end)
	     by (d (/ (- end begin) no-bins)))
	(setf (slot-value new-model p-name) x)
	(for int = (* 
		    (exp (funcall log-prior))
		    (integrate-over new-model data parameters-to-marginalize
				    :integration-function-creator integration-function-creator)))
	(sum int into integral)
	(collect (list x int) into vals)
	(finally
	    (return (if normalize
			(mapcar #'(lambda (x/y) (list (first x/y) (/ (second x/y) integral))) vals)
			vals)))))))
