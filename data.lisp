(in-package #:bayesian-analysis)

;;; classes

(defclass data ()
  ((no-data-points :initarg :no-data-points :accessor no-data-points)
   (independent-parameters :accessor independent-parameters :initarg :independent-parameters)
   (dependent-parameters :accessor dependent-parameters :initarg :dependent-parameters)
   (error-parameters :accessor error-parameters :initarg :error-parameters )
   (no-independent-parameters :accessor no-independent-parameters :initarg :no-independent-parameters)
   (no-dependent-parameters :accessor no-dependent-parameters :initarg :no-dependent-parameters)
   (descriptions :accessor descriptions :initarg :descriptions :initform '())
   (no-error-parameters :accessor no-error-parameters :initarg :no-error-parameters ))
  (:documentation "Super class for implementing types that represent
  data sets. DATA types' slots are used for keeping information about
  the subclass."))


(defclass parameter-description ()
  ((name :initarg :name :accessor name 
	 :initform (error "Must initialize name."))
   (textual-descriptoin :accessor textual-descriptoin :initarg :textual-descriptoin
			:initform "")))

;;; conditions


;;; api



(defun make-all-slots (independent-parameters dependent-parameters error-parameters)
  "Internal; create slot names and description for no, symbols, or
lists given to defmacro form."
  (labels ((mk-slot-names (symbol no-params-or-params)
	     (cond
	       ((and (integerp no-params-or-params) (>= no-params-or-params 0))
		(iter
		  (for i from 1 to no-params-or-params)
		  (collect (alexandria:symbolicate symbol (format nil "~D" i)) into name)
		  (collect "" into docs)
		  (finally (return (values name docs)))))
	       ((and (listp no-params-or-params)
		     (= 2 (length no-params-or-params))
		     (stringp (second no-params-or-params)))
		(values (list (first no-params-or-params))
			(list (second no-params-or-params))))
	       ((listp no-params-or-params)
		(iter
		  (for p in no-params-or-params)
		  (for name = (if (listp p) (first p) p))
		  (for doc = (if (and (listp p)
				      (stringp (second p)))
				 (second p)
				 ""))
		  (collect name into names)
		  (collect doc into docs)
		  (finally (return (values names docs)))))
	       ((symbolp no-params-or-params)
		(values (list no-params-or-params) (list "")))
	       (t (error "Parameter value must either be a number, symbol or list.")))))
    (let+ (((&values x-slots x-docs) (mk-slot-names 'x- independent-parameters))
	   ((&values y-slots y-docs) (mk-slot-names 'y- dependent-parameters))
	   ((&values error-slots error-docs) (mk-slot-names 's- error-parameters))
	   (all-slots (append x-slots y-slots error-slots))
	   (all-docs (append x-docs y-docs error-docs)))
      (values x-slots y-slots error-slots all-slots all-docs))))


(defun make-class-def (name all-slots)
  "Internal; create class definition based on parameters given to
defmacro form."
  `(defclass ,name (bayesian-analysis:data)
     ,(iter
	(for s in all-slots)
	(for key-x = (alexandria:make-keyword s))
	(collect `(,s :initarg ,key-x :accessor ,s)))))



(defun --init--check-congruent-data (object slot-names)
  "Upon initialization of a data object, this method will check if all
of the parameters satisfy the following conditions:

- all parameters must be represented by simple arrays with elements of
  type double-float.

- all parameter arrays must be of the same size."
  (let+ ((arrays (mapcar #'(lambda (s) (slot-value object s)) slot-names)))
    (map nil
	 #'(lambda (a s)
	     (if (not (and (typep a 'simple-array)
			   (eql 'double-float (array-element-type a))))
		 (error 'wrong-data-type
			:format-control
			"All parameters must be given as SIMPLE-ARRAYs
of type DOUBLE-FLOAT. Parameter ~a is of type: ~a"
			:format-arguments (list s (type-of a)))))
	 arrays slot-names))
  (setf (slot-value object 'no-data-points)
	(length (slot-value object (first slot-names)))))

(defun make-init-from-source-method (name init-from-source-body object-var-name
				     source-type source-var-name
				     all-slots)
  "Internal; create initialize-from-source defmethod for datatype
based on information given in defmacro form."
  (alexandria:with-gensyms (type)
    `(defmethod bayesian-analysis:initialize-from-source ((,type (eql ',name))
							  (,source-var-name ,source-type))
       (let-plus:let+ ((,object-var-name (make-instance ',name))
		       ((let-plus:&slots ,@all-slots) ,object-var-name))
	 (progn
	   ,@init-from-source-body)
	 (--init--check-congruent-data ,object-var-name ',all-slots)
	 ,object-var-name))))


(defun make-init-instance-method (name all-slots all-docs x-slots y-slots error-slots)
  "Internal; create initialize-instance defmethod for data based on
information given in defmacro form."
  `(defmethod initialize-instance :after ((object ,name) &key)
     (let-plus:let+ (((&slots no-independent-parameters independent-parameters
			      no-dependent-parameters dependent-parameters
			      no-error-parameters error-parameters
			      descriptions
			      ,@all-slots) object))
       (setf no-independent-parameters ,(length x-slots) no-dependent-parameters ,(length y-slots)
	     no-error-parameters ,(length error-slots) independent-parameters ',x-slots
	     dependent-parameters ',y-slots error-parameters ',error-slots
	     descriptions (mapcar #'(lambda (n d) (make-instance 'parameter-description
							    :name n
							    :textual-descriptoin d))
				  ',all-slots ',all-docs)))))

(defun make-get-parameter-names-method (type-name all independent dependent error)
  "Internal; create accesor defmethods for information about the
represented data type."
  (alexandria:with-gensyms (type)
    `((defmethod bayesian-analysis:get-all-data-slots ((,type (eql ',type-name)))
	',all)
      (defmethod bayesian-analysis:get-independent-parameters ((,type (eql ',type-name)))
	',independent)
      (defmethod bayesian-analysis:get-dependent-parameters ((,type (eql ',type-name)))
	',dependent)
      (defmethod bayesian-analysis:get-error-parameters ((,type (eql ',type-name)))
	',error))))

(defun macro-input-checks (y-slots error-slots)
  "Internal; check that the number of error-slots match the number of
dependent data slots."
  (cond
    ((and error-slots (not (= (length y-slots) (length error-slots))))
     (error 'wrong-number-of-arguments
	    :explanation "Need as many error parameters as dependent parameters."))))




(defmacro define-data-class (name
			     independent-parameters
			     dependent-parameters
			     error-parameters
			     (object-var-name (source-var-name source-type))
			     &body init-from-source-body)
  "Macro to define a data type to be used for likelihood calculations.
NAME is the type of the data class.  INDEPENDENT-PARAMETERS,
DEPENDENT-PARAMETERS, and ERROR-PARAMETERS can be given independently
as:

- a positive integer I, in which case the class slots representing
  independent, dependent, and/or error parameters respectively will be
  labeled X-1 .. X-I, Y-1 .. Y-I, S-1 .. S-I

- a single symbol that will be used for independent, dependent, and/or
  error parameters 

- a list containing a single symbol and a string: {DATUM,
  DESCRIPTION}, in which case DATUM will be used to name the parameter
  and DESCRIPTION.

The number of error parameters has to match the number of dependent
parameters. The descriptions, if given, will be stored and are
accessible at run-time.

Using parameters OBJECT-VAR-NAME, SOURCE-VAR-NAME, and SOURCE-TYPE,
the method bayesian-analysis:initialize-from-source will be
implemented using code given in INIT-FROM-SOURCE-BODY. All parameter
slots will be visible (as per with-slots) to the forms
INIT-FROM-SOURCE-BODY in addition to the object initialized which will
be bound to the variable named OBJECT-VAR-NAME. The method will be
specialized on (eql 'NAME) and an object of type SOURCE-TYPE, given to
the forms in INIT-FROM-SOURCE-BODY as method parameter
SOURCE-VAR-NAME.  INIT-FROM-SOURCE-BODY is expected to fill all
parameters with SIMPLE-ARRAYs of type DOUBLE-FLOAT (this is checked
automatically after INIT-FROM-SOURCE-BODY has been executed).
"
  (let+ (((&values x-slots y-slots error-slots all-slots all-descriptions)
	  (make-all-slots independent-parameters dependent-parameters error-parameters)))
    (macro-input-checks y-slots error-slots)
    `(eval-when (:compile-toplevel :load-toplevel)
       ,(make-class-def name all-slots)
       ,(make-init-instance-method name all-slots all-descriptions x-slots y-slots error-slots)
       ,(make-init-from-source-method name init-from-source-body object-var-name
				      source-type source-var-name all-slots)
       ,@(make-get-parameter-names-method name all-slots x-slots y-slots error-slots))))



