(in-package #:bayesian-analysis)



(defun model-function-name (model-name)
  (alexandria:symbolicate model-name '-model-function))

(defun make-model-function-access-functions (model-name model-parameters
					     function-independent-parameters
					     model-function-body data-data-type)
  (let+ ((model-function-name (model-function-name model-name))
	 (independent-parameters (get-independent-parameters data-data-type))
	 (dependent-parameters (get-dependent-parameters data-data-type))
	 (error-parameters (get-error-parameters data-data-type)) 
	 (f_i-name (alexandria:symbolicate 'f_i/ model-function-name))
	 (y_i-name (alexandria:symbolicate 'y_i/ model-function-name))
	 (err_i-name (alexandria:symbolicate 'err_i/ model-function-name))
	 (y_i-f_i-name (alexandria:symbolicate 'y_i-f_i/ model-function-name))
	 (sqrt/chi^2-2ln[p]-name (alexandria:symbolicate 'sqrt/chi^2-2ln[p]/ model-function-name))
	 (y_i-f_i/err_i-name (alexandria:symbolicate 'y_i-f_i/err_i/ model-function-name))
	 (model-function-code
	  `(defun ,model-function-name (,@function-independent-parameters model-object)
	     (declare (ignorable ,@function-independent-parameters)
		      (type double-float ,@function-independent-parameters)
		      (type ,model-name model-object))
	     (let-plus:let+ (((let-plus:&slots ,@model-parameters
					       ba:new-sample? ba:cache) model-object))
	       (declare (type double-float ,@model-parameters)
			(type (simple-array double-float) ba:cache))
	       (progn ,@model-function-body))))
	 (f_i-code
	  `(defun ,f_i-name (i model-object data-object)
	     (declare (type ,model-name model-object)
		      (type ,data-data-type data-object)
		      (type (integer 0) i))
	     (with-slots (,@independent-parameters) data-object
	       (declare (type (simple-array double-float) ,@independent-parameters))
	       (,model-function-name ,@(mapcar #'(lambda (p) `(aref ,p i )) independent-parameters) model-object))))
	 (y_i-code
	  `(defun ,y_i-name (i data-object)
	     ,@(if dependent-parameters
		   `((declare (type ,data-data-type data-object)
			      (type (integer 0) i))
		     (with-slots (,@dependent-parameters) data-object
		       (declare (type (simple-array double-float) ,@dependent-parameters))
		       (aref ,@dependent-parameters i)))
		   `((declare (ignore i data-object))
		     0d0))))
	 (err_i-code
	  `(defun ,err_i-name (i data-object)
	     ,@(if error-parameters
		   `((declare (type ,data-data-type data-object)
			      (type (integer 0) i))
		     (with-slots (,@error-parameters) data-object
		       (declare (type (simple-array double-float) ,@error-parameters))
		       (aref ,@error-parameters i)))
		   `((declare (ignore i data-object)) 0d0))))
	 (y_i-f_i-code
	  `(defun ,y_i-f_i-name (i model-object data-object)
	     ,@(if dependent-parameters
		   `((declare (type ,model-name model-object)
			      (type ,data-data-type data-object)
			      (type (integer 0) i))
		     (with-slots (,@independent-parameters ,@dependent-parameters) data-object
		       (declare (type (simple-array double-float)
				      ,@independent-parameters ,@dependent-parameters))
		       (- (aref ,@dependent-parameters i)
			  (,model-function-name ,@(mapcar #'(lambda (p) `(aref ,p i )) independent-parameters)
						model-object))))
		   `((declare (ignore i model-object data-object)) 0d0))))
	 ;; (sqrt/chi^2-2ln[p]-code
	 ;;  `(defun ,sqrt/chi^2-2ln[p]-name (i model-object data-object)
	 ;;    ,@(if dependent-parameters
	 ;; 	   `((declare (type ,model-name model-object)
	 ;; 		      (type ,data-data-type data-object)
	 ;; 		      (type (integer 0) i))
	 ;; 	     (with-slots (,@independent-parameters ,@dependent-parameters) data-object
	 ;; 	       (declare (type (simple-array double-float)
	 ;; 			      ,@independent-parameters ,@dependent-parameters))
	 ;; 	       (let-plus:let+ (((&slots log-of-all-priors) ,model-object)
	 ;; 			       (Q (- (aref ,@dependent-parameters i)
	 ;; 				     (,model-function-name
	 ;; 				      ,@(mapcar #'(lambda (p) `(aref ,p i )) independent-parameters)
	 ;; 				      model-object))))
	 ;; 		 (sqrt (- (* Q Q)
	 ;; 			  (funcall log-of-all-priors))))))
	 ;; 	   `((declare (ignore i model-object data-object)) 0d0))))
	 (y_i-f_i/err_i-code
	  `(defun ,y_i-f_i/err_i-name (i model-object data-object)
	     ,@(if (and dependent-parameters error-parameters)
		   `((declare (type ,model-name model-object)
			      (type ,data-data-type data-object)
			      (type (integer 0) i))
		     (with-slots (,@independent-parameters ,@dependent-parameters ,@error-parameters)
			 data-object
		       (declare (type (simple-array double-float)
				      ,@independent-parameters ,@dependent-parameters ,@error-parameters))
		       (/ (- (aref ,@dependent-parameters i)
			     (,model-function-name ,@(mapcar #'(lambda (p) `(aref ,p i )) independent-parameters)
						   model-object))
			  (aref ,@error-parameters i))))
		   `((declare (ignore i model-object data-object)) 0d0)))))
    (unless (>= (length independent-parameters) 1)
      (error "Need at least one dependent parameter."))
    (values model-function-code
	    f_i-code y_i-code err_i-code
	    y_i-f_i-code y_i-f_i/err_i-code
	    f_i-name y_i-name err_i-name
	    y_i-f_i-name y_i-f_i/err_i-name)))
