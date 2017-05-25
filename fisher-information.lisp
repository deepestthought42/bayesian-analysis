(in-package #:bayesian-analysis)


(defclass fisher-information ()
  ((matrix :initarg :matrix :accessor matrix 
	   :initform (error "Must initialize matrix."))
   (parameter-names :initarg :parameter-names :accessor parameter-names 
		    :initform (error "Must initialize parameter-names."))
   (no-parameters :initarg :no-parameters :accessor no-parameters 
		  :initform (error "Must initialize no-parameters."))
   (weighted-by-errors :initarg :weighted-by-errors :accessor weighted-by-errors 
		       :initform (error "Must initialize weighted-by-errors."))
   (determinant :initarg :determinant :accessor determinant 
		:initform (error "Must initialize determinant."))))


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
			 (min
			  1d-3
			  (* (expt epsilon 0.25d0)
			     (slot-value model param)))))))))


(defun hessian (func model params.delta &optional (sign 1d0))
  "Calculate the negative hessian matrix for FUNC, where FUNC is a
function object (closure) that depends on MODEL. PARAMS.DELTA is a
list of (PARAMETER-SLOT DELTA), where PARAMETER-SLOT is the name of a
slot that was marginalized and DELTA is the optimal delta for that
variable.

The calculation is based on h_{3,j,k}, eq. 9, in: 


Statistical Applications of the Complex-Step Method of Numerical Differentiation
Author(s): Martin S. Ridout
Source: The American Statistician, Vol. 63, No. 1 (Feb., 2009), pp. 66-74

For reference, the formula itself is:

\begin{equation}
\label{eq:approximate-hessian}
h_{j,k}=-\frac{1}{4\delta_j\delta_k}
        \left\{\left[
                 f \left(\mathbf{\theta}+\delta_{j}\mathbf{e}_j + \delta_k\mathbf{e}_k \right)
                 - f\left(\mathbf{\theta}+\delta_{j}\mathbf{e}_j - \delta_k\mathbf{e}_k \right)
                \right] 
                -
                \left[
                 f \left(\mathbf{\theta}-\delta_{j}\mathbf{e}_j + \delta_k\mathbf{e}_k \right)
                 - f\left(\mathbf{\theta}-\delta_{j}\mathbf{e}_j - \delta_k\mathbf{e}_k \right)
                \right] 
        \right\}
\end{equation}

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
		 (d param-j delta-j)
		 (d param-k delta-k)
		 (setf a (funcall func)) ; f(t + d_je_j + d_ke_k)
		 (d param-k (- (* 2d0 delta-k)))
		 (setf b (funcall func)) ; f(t + d_je_j - d_ke_k)
		 (d param-j (- (* 2d0 delta-j)))
		 (setf d (funcall func)) ; f(t - d_je_j - d_ke_k)
		 (d param-k (* 2d0 delta-k))
		 (setf c (funcall func)) ; f(t - d_je_j + d_ke_k)
		 (* sign
		    (/ (- (- a b)
			  (- c d))
		       (* 4d0 delta-j delta-k))))))
      (iter
	(for j from 0 below dim)
	(iter
	  (for k from j below dim)
	  (let+ ((grad (h-j-k (param j) (delta j)
			      (param k) (delta k))))
	    (setf (aref ret-val j k) grad
		  (aref ret-val k j) grad))))
      ret-val)))


(defmethod get-fisher-information-matrix ((result nlopt-result) &key)
  "See generic function for documentation."
  (let+ (((&slots model data likelihood) result)
	 ((&slots log-of-all-priors) model)
	 ((&slots varying/log-of-likelihood constant/log-of-likelihood) likelihood))
    (labels ((fun ()
	       (+ (funcall varying/log-of-likelihood)
		  (funcall constant/log-of-likelihood)
		  (funcall log-of-all-priors))))
      (hessian #'fun model (get-optimal-delta model) -1d0))))


#+nil
(progn
  (defun test-fisher-information-matrix ()
    (let+ ((model (make-instance 'linear))
	   (data (initialize-from-source '1d-data t))
	   (result
	    (find-optimum (make-instance 'levenberg-marquardt) model data))
	   (matrix (get-fisher-information-matrix result)))
      matrix))

  (test-fisher-information-matrix))


