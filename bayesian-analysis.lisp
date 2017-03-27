;;;; bayesian-analysis.lisp

(in-package #:bayesian-analysis)

;;; "bayesian-analysis" goes here. Hacks and glory await!



;;; api

(defgeneric get-prior-for-parameter (model parameter)
  (:documentation "Return the priors for parameter PAREMETER 
=> log-of-prior, varying/log-of-prior, constant/log-of-prior"))

;;; data api

(defgeneric initialize-from-source (data-type source)
  (:documentation "Initialize an object of type DATA-TYPE form soure SOURCE."))

(defgeneric get-all-data-slots (data-type)
  (:documentation "Get all slots known for type DATA-TYPE."))

(defgeneric get-independent-parameters (type)
  (:documentation "Get all independent parameters for type DATA-TYPE."))

(defgeneric get-dependent-parameters (type)
  (:documentation "Get all dependent parameters for type DATA-TYPE."))

(defgeneric get-error-parameters (type)
  (:documentation "Get all error parameters for type DATA-TYPE."))

(defgeneric get-parameter-description-text (object parameter)
  (:documentation "Get the textual description of the parameter
  PARAMETER for data object OBJECT."))


;;; conditions

(define-condition unknown-prior-type (error)
  ((type-not-known :accessor type-not-known :initarg :type-not-known
		   :initform :unspecified)))

(define-condition initial-parameter-out-of-range ()
  ((parameter :accessor parameter :initarg :parameter :initform :unknown)))

(define-condition parameter-out-of-range (error)
  ((parameter :accessor parameter :initarg :parameter :initform :unknown)))

(define-condition unknown-sampling-type (error)
  ((type-not-known :accessor type-not-known :initarg :type-not-known
		   :initform :unspecified)))


(define-condition incongruent-data (simple-condition) ())
(define-condition wrong-data-type (simple-condition) ())
(define-condition unknown-parameter (simple-condition) ())
(define-condition wrong-number-of-arguments (program-error)
  ((explanation :accessor explanation :initarg :explanation :initform "")
   (offending-symbol :accessor offending-symbol :initarg :offending-symbol :initform nil)))


(define-condition no-parameters-to-marginalize (simple-condition) ())
