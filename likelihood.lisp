(in-package #:bayesian-analysis)



(defclass likelihood ()
  ((constant/log-of-likelihood :accessor constant/log-of-likelihood :initarg
			       :constant/log-of-likelihood
			       :initform (constantly 0d0))
   (varying/log-of-likelihood :accessor varying/log-of-likelihood :initarg :varying/log-of-likelihood
			      :initform (constantly 0d0))))

(defgeneric likelihood (likelihood)
  (:method ((l likelihood))
    (exp (+ (funcall (constant/log-of-likelihood l))
	    (funcall (varying/log-of-likelihood l))))))


(defparameter *likelihood-types*
  '((:d_i=f_i+gaussian_error_i_unequal_sigma)
    (:d_i=f_i+gaussian_error_1_equal_sigma)))



(defun %check-likelihood-params (likelihood-type equal-sigma-parameter)
  (case likelihood-type
    (:d_i=f_i+gaussian_error_i_unequal_sigma)
    (:d_i=f_i+gaussian_error_1_equal_sigma
     (if (not equal-sigma-parameter)
	 (error "Need to provide parameter definition for EQUAL-SIGMA-PARAMETER, 
when using :d_i=f_i+gaussian_error_1_equal_sigma type likelihood.")))))


(defun create-likelihood-functions/gaussian/i-known-errors (model-function-name
							    model-object
							    data-object)
  (let+ ((xs (slot-value data-object (first (independent-parameters data-object))))
	 (ys (slot-value data-object (first (dependent-parameters data-object))))
	 (errors (slot-value data-object (first (error-parameters data-object))))
	 (no-data-points (no-data-points data-object)))
    (labels ((varying ()
	       (iter
		 (with Q = 0d0)
		 (for i from 0 below no-data-points)
		 (for x = (aref xs i))
		 (for y = (aref ys i))
		 (for err = (aref errors i))
		 (for d_i-y_i = (- y (funcall model-function-name x model-object)))
		 (for sqrt-q = (/ d_i-y_i err))
		 (incf Q (* sqrt-q sqrt-q))
		 (finally (return (- (/ Q 2d0))))))
	     (constant ()
	       (iter
		 (with retval = 1d0)
		 (for i from 0 below no-data-points)
		 (setf retval (* retval (aref errors i)))
		 (finally
		  (return
		    (log (/ 1 (* (expt (* 2 pi) (/ no-data-points 2)) retval))))))))
      (make-instance 'likelihood :varying/log-of-likelihood #'varying
				 :constant/log-of-likelihood #'constant))))

(defun create-likelihood-functions/gaussian/1-unknown-error (model-function-name
							     model-object
							     data-object
							     equal-sigma-parameter
							     amplitude-parameter)
  (let+ ((xs (slot-value data-object (first (independent-parameters data-object))))
	 (ys (slot-value data-object (first (dependent-parameters data-object))))
	 (no-data-points (no-data-points data-object))
	 (N (coerce no-data-points 'double-float)))
    (labels ((varying ()
	       (iter
		 (with Q = 0d0)
		 (for i from 0 below no-data-points)
		 (for x = (aref xs i))
		 (for y = (aref ys i))
		 (for d_i-y_i = (- y (funcall model-function-name x model-object)))
		 (incf Q (* d_i-y_i d_i-y_i))
		 ;; if this turns out to be a performance problem, we
		 ;; can put it into a base class
		 (finally
		  (let ((sigma (slot-value model-object equal-sigma-parameter))
			(amplitude (slot-value model-object amplitude-parameter)))
		    (return (- (+ (* N (log (* sigma amplitude)))
				  (/ Q (* 2d0 sigma sigma)))))))))
	     (constant ()
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




(defun make-likelihood-initializer (model-name likelihood-type
				    model-function-name
				    data-dependent-parameters
				    data-error-parameters
				    data-type
				    equal-sigma-parameter
				    amplitude-parameter)
  (alexandria:with-gensyms (data-object-name model-object-name)
    `(defmethod bayesian-analysis:initialize-likelihood ((,model-object-name ,model-name)
							 (,data-object-name ,data-type))
       ,(case likelihood-type
	  (:d_i=f_i+gaussian_error_i_unequal_sigma
	   (if (not (= 1 (length data-dependent-parameters) (length data-error-parameters)))
	       (error "For a simple d_i = f_i + e_i, only one error
	 parameter and one dependent parameter is supported."))
	   `(create-likelihood-functions/gaussian/i-known-errors #',model-function-name
								 ,model-object-name
								 ,data-object-name))
	  (:d_i=f_i+gaussian_error_1_equal_sigma
	   (if (not (and (= 1 (length data-dependent-parameters))
			 (=  (length data-error-parameters))))
	       (error "For a simple d_i = f_i + e_i, only one error
	 parameter and one dependent parameter is supported."))
	   `(create-likelihood-functions/gaussian/1-unknown-error #',model-function-name
								  ,model-object-name
								  ,data-object-name
								  ',equal-sigma-parameter
								  ',amplitude-parameter))))))

 
