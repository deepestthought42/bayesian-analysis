(in-package #:bayesian-analysis)


(defclass fisher-information ()
  ((covariance-matrix-from-gsl :initarg :covariance-matrix-from-gsl :accessor covariance-matrix-from-gsl 
			       :initform (error "Must initialize covariance-matrix-from-gsl."))
   (matrix :initarg :matrix :accessor matrix 
	   :initform (error "Must initialize matrix."))
   (parameter-names :initarg :parameter-names :accessor parameter-names 
		    :initform (error "Must initialize parameter-names."))
   (no-parameters :initarg :no-parameters :accessor no-parameters 
		  :initform (error "Must initialize no-parameters."))
   (weighted-by-errors :initarg :weighted-by-errors :accessor weighted-by-errors 
		       :initform (error "Must initialize weighted-by-errors."))
   (determinant :initarg :determinant :accessor determinant 
		:initform (error "Must initialize determinant."))))

(defgeneric get-fisher-information-matrix (result &key))

;; fixmee: need to reenter this in asd


(defun get-optimal-delta (model &optional (epsilon long-float-epsilon epsilon-given-p))
  (let+ (((&slots model-parameters-to-marginalize) model))
    (iter
      (for param in model-parameters-to-marginalize)
      ;; fixme: should look up what happens if the value is below the
      ;; machine accuracy
      (collect (list param
		     (if epsilon-given-p
			 epsilon
			 (* (expt epsilon 0.25d0)
			    (slot-value model param))))))))


(defun negative-hessian (func model params.delta)
  "Calculate the negative hessian matrix for FUNC, where FUNC is a
function object (closure) that depends on MODEL. PARAMS.DELTA is a
list of (PARAMETER-SLOT DELTA), where PARAMETER-SLOT is the name of a
slot that was marginalized and DELTA is the optimal delta for that
variable.
"
  (let+ ((dim (length params.delta))
	 ;; fixmee: type information here
	 (ret-val (make-array (list dim dim))))
    (labels ((param (i) (first (nth i params.delta)))
	     (delta (i) (second (nth i params.delta)))
	     (d (param delta)
	       (incf (slot-value model param) delta))
	     (h-j-k (param-j delta-j param-k delta-k)
	       (let ((a 0d0) (b 0d0)
		     (c 0d0) (d 0d0))
		 (d param-j delta-j) (d param-k delta-k)
		 (setf a (funcall func))
		 (d param-k (- (* 2d0 delta-k)))
		 (setf b (funcall func))
		 (d param-j (- (* 2d0 delta-j)))
		 (setf d (funcall func))
		 (d param-k (* 2d0 delta-k))
		 (setf c (funcall func))
		 (/ (- (- a b)
		       (- c d))
		    (* 4d0 delta-j delta-k)))))
      (iter
	(for j from 0 below dim)
	(iter
	  (for k from j below dim)
	  (let+ ((grad (h-j-k (param j) (delta j)
			      (param k) (delta k))))
	    (setf (aref ret-val j k) grad
		  (aref ret-val k j) grad))))
      ret-val)))



(defmethod get-fisher-information-matrix ((result levenberg-marquardt-parameter-result)
					  &key (count-class 1))
  ;; fixme: not taking count classes into account
  (declare (ignore count-class))
  (let+ (((individual-result &rest &ign) (individual-fit-results result))
	 ((&slots fit-model weighted-by-errors) result)
	 ((&slots gf:no-fit-params gf:slots-to-fit) fit-model)
	 (covariance-matrix (gf:result-covar individual-result))
	 ((&values fisher-matrix determinant) (math-utils:invert-matrix covariance-matrix)))
    (make-instance 'fisher-information
		   :covariance-matrix-from-gsl covariance-matrix
		   :no-parameters gf:no-fit-params
		   :parameter-names gf:slots-to-fit
		   :matrix fisher-matrix
		   :weighted-by-errors weighted-by-errors
		   :determinant determinant)))

