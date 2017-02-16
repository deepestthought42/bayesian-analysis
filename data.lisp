(in-package #:bayesian-analysis)


(defclass data ()
  ((no-data-points :initarg :no-data-points :accessor no-data-points)
   (independent-parameters :accessor independent-parameters :initarg :independent-parameters)
   (dependent-parameters :accessor dependent-parameters :initarg :dependent-parameters)
   (error-parameters :accessor error-parameters :initarg :error-parameters )
   (no-independent-parameters :accessor no-independent-parameters :initarg :no-independent-parameters)
   (no-dependent-parameters :accessor no-dependent-parameters :initarg :no-dependent-parameters)
   (no-error-parameters :accessor no-error-parameters :initarg :no-error-parameters )))



(defun make-all-slots (independent-parameters dependent-parameters error-parameters)
  (labels ((mk-slot-names (symbol no-params-or-params)
	     (cond
	       ((numberp no-params-or-params)
		(iter
		  (for i from 1 to no-params-or-params)
		  (collect (alexandria:symbolicate symbol (format nil "~D" i)))))
	       ((listp no-params-or-params) no-params-or-params)
	       ((symbolp no-params-or-params) (list no-params-or-params))
	       (t (error "Parameter value must either be a number, symbol or list.")))))
    (let+ ((x-slots (mk-slot-names 'x- independent-parameters))
	   (y-slots (mk-slot-names 'y- dependent-parameters))
	   (error-slots (mk-slot-names 's- error-parameters))
	   (all-slots (append x-slots y-slots error-slots)))
      (values x-slots y-slots error-slots all-slots))))


(defun make-class-def (name all-slots)
  `(defclass ,name (bayesian-analysis:data)
     ,(iter
	(for s in all-slots)
	(for key-x = (alexandria:make-keyword s))
	(collect `(,s :initarg ,key-x :accessor ,s)))))

(defun make-init-from-source-method (name init-from-source-body object-var-name
				     source-type source-var-name
				     source-var-name-f-data-points get-no-data-points-body
				     x-slots y-slots error-slots all-slots)
  `(defmethod bayesian-analysis:initialize-from-source ((eql ,name)
							(,source-var-name ,source-type))
     (let-plus:let+ ((,object-var-name (make-instance ',name))
		     ((let-plus:&slots ,@all-slots
				       no-independent-parameters independent-parameters
				       no-dependent-parameters dependent-parameters
				       no-error-parameters error-parameters
				       no-data-points) ,object-var-name))
       (setf no-independent-parameters ,(length x-slots) no-dependent-parameters ,(length y-slots)
	     no-error-parameters ,(length error-slots) independent-parameters ',x-slots
	     dependent-parameters ',y-slots error-parameters ',error-slots
	     no-data-points
	     (funcall #'(lambda (,source-var-name-f-data-points)
			  (declare (ignorable ,source-var-name-f-data-points))
			  (progn
			    ,@get-no-data-points-body))
		      ,source-var-name))
       (progn
	 ,@init-from-source-body))))

(defgeneric initialize-from-source (data-object source))

(defmacro define-data-class (name
			     independent-parameters
			     dependent-parameters
			     error-parameters
			     ((source-var-name-f-data-points) &body get-no-data-points-body)
			     (object-var-name (source-var-name source-type))
			     &body init-from-source-body)
  (let+ (((&values x-slots y-slots error-slots all-slots)
	  (make-all-slots independent-parameters dependent-parameters error-parameters)))
    `(progn
       ,(make-class-def name all-slots)
       ,(make-init-from-source-method name init-from-source-body object-var-name
				      source-type source-var-name
				      source-var-name-f-data-points get-no-data-points-body
				      x-slots y-slots error-slots all-slots))))





(define-data-class 1d-data x y s
    ((source) 0)
    (object (source t)))
