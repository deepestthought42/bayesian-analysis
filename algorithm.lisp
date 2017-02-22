(in-package #:bayesian-analysis)


(defclass algorithm () ())

(defclass mcmc-algorithm (algorithm)
  ((no-iterations :accessor no-iterations
		  :initarg :no-iterations :initform 1000)))


(defclass mcmc-parameter-result ()
  ((algorithm :initarg :algorithm :accessor algorithm 
	      :initform (error "Must initialize algorithm."))
   (no-accepted-iterations :initarg :no-accepted-iterations :accessor no-accepted-iterations 
			   :initform 0)
   (no-iterations :initarg :no-iterations :accessor no-iterations :initform 0)
   (data :initarg :data :accessor data 
	 :initform (error "Must initialize data."))
   (input-model :initarg :input-model :accessor input-model 
		:initform (error "Must initialize input-model."))
   (result-model :initarg :result-model :accessor result-model 
		 :initform (error "Must initialize result-model."))
   (iteration-accumulator :initarg :iteration-accumulator :accessor iteration-accumulator 
			  :initform (error "Must initialize iteration-accumulator."))))


(defgeneric solve-for-parameters (algorithm model data &key))


(defclass mcmc-accumulator ()
  ((parameter-array :initarg :parameter-array :accessor parameter-array 
		    :initform (error "Must initialize parameter-array."))
   (no-iterations :initarg :no-iterations :accessor no-iterations 
		  :initform (error "Must initialize no-iterations."))
   (marginalized-parameters :initarg :marginalized-parameters :accessor marginalized-parameters 
			    :initform (error "Must initialize marginalized-parameters."))
   (no-marginalized-parameters :accessor no-marginalized-parameters :initarg :no-marginalized-parameters)
   (save-current-parameters :initarg :save-current-parameters :accessor save-current-parameters 
			    :initform (error "Must initialize save-current-parameters."))
   (save-last-parameters :initarg :save-last-parameters :accessor save-last-parameters 
			 :initform (error "Must initialize save-last-parameters."))))

(defmethod initialize-instance :after ((object mcmc-accumulator) &key)
  (let+ (((&slots no-marginalized-parameters marginalized-parameters) object))
    (setf no-marginalized-parameters (length marginalized-parameters))))

(defgeneric initialize-accumulator (model no-iterations))
