(in-package #:bayesian-analysis)

(declaim (optimize (debug 3) (space 0) (safety 1) (speed 3)))


(defclass likelihood ()
  ((constant/log-of-likelihood :accessor constant/log-of-likelihood :initarg
			       :constant/log-of-likelihood
			       :initform (constantly 0d0))
   (varying/log-of-likelihood :accessor varying/log-of-likelihood :initarg :varying/log-of-likelihood
			      :initform (constantly 0d0))))

(defgeneric likelihood (likelihood)
  (:method ((l likelihood))
    #+sbcl (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
    (exp (+ (funcall (constant/log-of-likelihood l))
	    (funcall (varying/log-of-likelihood l))))))


(defparameter *likelihood-types*
  '((:d_i=f_i+gaussian_error_i_unequal_sigma)
    (:d_i=f_i+gaussian_error_1_equal_sigma)
    (:p_of_x_i=f_i_of_x_i)))



(defun %check-likelihood-params (likelihood-type equal-sigma-parameter)
  (case likelihood-type
    (:d_i=f_i+gaussian_error_i_unequal_sigma)
    (:p_of_x_i=f_i_of_x_i)
    (:d_i=f_i+gaussian_error_1_equal_sigma
     (if (not equal-sigma-parameter)
	 (error "Need to provide parameter definition for EQUAL-SIGMA-PARAMETER, 
when using :d_i=f_i+gaussian_error_1_equal_sigma type likelihood.")))))


(defun create-likelihood-functions/gaussian/i-known-errors (model-object data-object
							    y_i-f_i/err_i err_i)
  (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
  (let+ ((no-data-points (no-data-points data-object))
	 (N (coerce no-data-points 'double-float)))
    (declare (type (integer 0) no-data-points))
    (labels ((varying ()
	       (iter
		 (declare (type fixnum i no-data-points)
			  (type (double-float 0d0) Q sqrt-q)
			  (type (function (fixnum t t) double-float) y_i-f_i/err_i))
		 (with Q = 0d0)
		 (for i from 0 below no-data-points)
		 (for sqrt-q = (funcall y_i-f_i/err_i i model-object data-object))
		 (incf Q (* sqrt-q sqrt-q))
		 (finally (return (- (/ Q 2d0))))))
	     (constant ()
	       (iter
		 (declare (type double-float retval)
			  (type (double-float 0d0) N)
			  (type fixnum i no-data-points)
			  (type (function (fixnum t) double-float) err_i))
		 (with retval = 1d0)
		 (for i from 0 below no-data-points)
		 (setf retval (* retval (funcall err_i i data-object)))
		 (finally
		  (return
		    (log (the (double-float 0d0)
			      (/ 1 (* (expt (* 2d0 pi) (/ N 2d0)) retval)))))))))
      (declare (ftype (function () double-float) varying constant))
      (make-instance 'likelihood :varying/log-of-likelihood #'varying
				 :constant/log-of-likelihood #'constant))))

(defun create-likelihood-functions/gaussian/1-unknown-error (model-object data-object
							     equal-sigma-parameter y_i-f_i)
  (let+ ((no-data-points (no-data-points data-object))
	 (N (coerce no-data-points 'double-float)))
    (labels ((varying ()
	       (iter
		 (declare (type fixnum i no-data-points)
			  (type (double-float 0d0) Q d_i-y_i)
			  (type (function (fixnum t t) double-float) y_i-f_i))
		 (with Q = 0d0)
		 (for i from 0 below no-data-points)
		 (for d_i-y_i = (funcall y_i-f_i i model-object data-object))
		 (incf Q (* d_i-y_i d_i-y_i))
		 ;; if this turns out to be a performance problem, we
		 ;; can put it into a base class
		 (finally
		  (let ((sigma (slot-value model-object equal-sigma-parameter)))
		    (declare (type (double-float 0d0) sigma))
		    (return (- (+ (* N (log sigma)) (/ Q (* 2d0 sigma sigma)))))))))
	     (constant ()
	       (declare (type (double-float 0d0) N))
	       (expt (* 2d0 pi) (/ N 2d0))))
      (make-instance 'likelihood
		     :varying/log-of-likelihood
		     (if *debug-function*
			 #'(lambda ()
			     (let ((val (varying)))
			       (debug-out :info :likelihood
					  "Caluclated varying likelihood to be: ~f" val)
			       val))
			 #'varying)
		     :constant/log-of-likelihood #'constant))))

(defun create-likelihood-functions/direct-distribution (model-object
							data-object
							f_i)
  (let+ ((no-data-points (no-data-points data-object)))
    (labels ((varying ()
	       (iter
		 (declare (type fixnum i no-data-points)
			  (type double-float Q)
			  (type (function (fixnum t t) double-float) f_i))
		 (with Q = 0d0)
		 (for i from 0 below no-data-points)
		 (incf Q (log (funcall f_i i model-object data-object)))
		 (finally
		  (return Q))))
	     (constant () 0d0))
      (make-instance 'likelihood
		     :varying/log-of-likelihood
		     (if *debug-function*
			 #'(lambda ()
			     (let ((val (varying)))
			       (debug-out :info :likelihood
					  "Caluclated varying likelihood to be: ~f" val)
			       val))
			 #'varying)
		     :constant/log-of-likelihood #'constant))))




(defun make-likelihood-initializer (model-name likelihood-type
				    data-type
				    equal-sigma-parameter
				    f_i-name y_i-name err_i-name
				    y_i-f_i-name y_i-f_i/err_i-name)
  (declare (ignore y_i-name))
  (alexandria:with-gensyms (data-object-name model-object-name)
    `(defmethod bayesian-analysis:initialize-likelihood ((,model-object-name ,model-name)
							 (,data-object-name ,data-type))
       ,(let ()
	  #+sbcl (declare (sb-ext:muffle-conditions sb-ext:compiler-note))
	  (case likelihood-type
	    (:d_i=f_i+gaussian_error_i_unequal_sigma
	     `(create-likelihood-functions/gaussian/i-known-errors ,model-object-name
								   ,data-object-name
								   #',y_i-f_i/err_i-name
								   #',err_i-name))
	    (:d_i=f_i+gaussian_error_1_equal_sigma
	     `(create-likelihood-functions/gaussian/1-unknown-error ,model-object-name
								    ,data-object-name
								    ',equal-sigma-parameter
								    #',y_i-f_i-name))
	    (:p_of_x_i=f_i_of_x_i
	     `(create-likelihood-functions/direct-distribution ,model-object-name
							       ,data-object-name
							       #',f_i-name))
	    (t (error "Unknown type of likelihood: ~a" likelihood-type)))))))

 
