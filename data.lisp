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
  (alexandria:with-gensyms (type)
    `(defmethod bayesian-analysis:initialize-from-source ((,type (eql ',name))
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
	   ,@init-from-source-body)
	 (if (not (apply #'= (mapcar #'length (list ,@all-slots))))
	     (error "All data arrays need to be of the same size!"))
	 (setf no-data-points (length ,(first all-slots)))
	 ,object-var-name))))

(defun make-get-parameter-names-method (type-name all independent dependent error)
  (alexandria:with-gensyms (type)
    `((defmethod bayesian-analysis:get-all-data-slots ((,type (eql ',type-name)))
	',all)
      (defmethod bayesian-analysis:get-independent-parameters ((,type (eql ',type-name)))
	',independent)
      (defmethod bayesian-analysis:get-dependent-parameters ((,type (eql ',type-name)))
	',dependent)
      (defmethod bayesian-analysis:get-error-parameters ((,type (eql ',type-name)))
	',error))))

(defgeneric initialize-from-source (data-object source))
(defgeneric get-all-data-slots (type))
(defgeneric get-independent-parameters (type))
(defgeneric get-dependent-parameters (type))
(defgeneric get-error-parameters (type))

(defun input-checks (y-slots error-slots)
  (cond
    ((if (and error-slots (not (= (length y-slots) (length error-slots))))
	 (error "Need as many error parameters as dependent parameters.")))))

(defmacro define-data-class (name
			     independent-parameters
			     dependent-parameters
			     error-parameters
			     ((source-var-name-f-data-points) &body get-no-data-points-body)
			     (object-var-name (source-var-name source-type))
			     &body init-from-source-body)
  (let+ (((&values x-slots y-slots error-slots all-slots)
	  (make-all-slots independent-parameters dependent-parameters error-parameters)))
    (input-checks y-slots error-slots)
    `(progn
       ,(make-class-def name all-slots)
       ,(make-init-from-source-method name init-from-source-body object-var-name
				      source-type source-var-name
				      source-var-name-f-data-points get-no-data-points-body
				      x-slots y-slots error-slots all-slots)
       ,@(make-get-parameter-names-method name all-slots x-slots y-slots error-slots))))







