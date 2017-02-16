(in-package #:bayesian-analysis)


(defclass data ()
  ((y :initarg :y :accessor y 
      :initform (error "Must initialize y."))
   (no-data-points :initarg :no-data-points :accessor no-data-points)))


(defclass data-with-sigma (data)
  ((sigma :initarg :sigma :accessor sigma 
	  :initform (error "Must initialize sigma."))))



(defgeneric initialize-from-source (data-object source))

(defmacro define-data-class (name no-params with-sigma
			     (object-var-name (source-var-name source-type))
			     &body init-from-source-body)
  (let+ ((base-slots (if with-sigma '(ba:y ba:sigma) '(ba:y)))
	 (x-slots (iter
		    (for i from 1 to no-params)
		    (collect (alexandria:symbolicate 'x- (format nil "~D" i)))))
	 (all-slots (append base-slots x-slots)))
    `(progn
       (defclass ,name (,(if with-sigma
			     'bayesian-analysis:data-with-sigma
			     'bayesian-analysis:data))
	 ,(iter
	    (for x in x-slots)
	    (for key-x = (alexandria:make-keyword x))
	    (collect `(,x :initarg ,key-x :accessor ,x
			  :initform (error ,(format nil "Must initialize ~a." x))))))
       (defmethod bayesian-analysis:initialize-from-source ((,object-var-name ,name)
							    (,source-var-name ,source-type))
	 (let-plus:let+ (((let-plus:&slots ,@all-slots) ,object-var-name))
	   (progn
	   ,@init-from-source-body))))))










