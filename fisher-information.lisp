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

(defmethod get-fisher-information-matrix ((result result) &key (count-class 1))
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
