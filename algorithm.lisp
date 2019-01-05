(in-package #:bayesian-analysis)


(defclass algorithm () ())


(defclass optimization-result ()
  ((algorithm :initarg :algorithm :accessor algorithm 
	      :initform (error "Must initialize algorithm."))
   (data :initarg :data :accessor data 
	 :initform (error "Must initialize data."))
   (input-model :initarg :input-model :accessor input-model 
		:initform (error "Must initialize input-model."))))





(defgeneric find-optimum (algorithm model data &key))



(defgeneric initialize-accumulator (model no-iterations))


(defgeneric bin-parameter-values (result parameter &key start end))
(defgeneric get-parameter-results (result &key start end))
(defgeneric get-optimal-values (result &key start end))


(defun %calculate-confidance (array normalized-to confidance-level)
  (let+ ((sorted (sort (copy-seq array) #'> :key #'second)))
    (iter
      (for (x counts) in-sequence sorted)
      (sum counts into summed)
      (maximize x into max)
      (minimize x into min)
      (if (>= summed (* normalized-to confidance-level))
	  (return (values min max)))
      (finally (error "Huh ? Couldn't calculate confidance interval.")))))



(defun %get-bin-width/no-bins (parameter-array pos start end no-bins
			       &key round)
  (let+ (((&values min max unrounded)
	  (iter
	    (for i from start below end)
	    (for val = (aref parameter-array pos i))
	    (maximize val into max)
	    (minimize val into min)
	    (finally
	         (if (<= max min)
		     (error "Couldn't calculate bin width."))
		 (return (values min max (/ (- max min) no-bins))))))
	 (bin-width (if round
			(expt 10 (floor (log unrounded 10)))
			unrounded))
	 (no-bins (ceiling (/ (- max min) bin-width))))
    (values (coerce bin-width 'double-float) no-bins)))




(defclass optimized-parameters ()
  ((model :initarg :model :accessor model 
	  :initform (error "Must initialize model."))
   (parameter-infos :initarg :parameter-infos :accessor parameter-infos 
		    :initform (error "Must initialize parameter-infos."))
   (algorithm-result :initarg :algorithm-result :accessor algorithm-result 
		     :initform (error "Must initialize algorithm-result."))
   (data :initarg :data :accessor data 
	 :initform (error "Must initialize data."))))




(defclass parameter-distribution ()
  ((name :initarg :name :accessor name 
	 :initform (error "Must initialize name."))
   (median :initarg :median :accessor median 
	   :initform (error "Must initialize median."))
   (confidence-level :initarg :confidence-level :accessor confidence-level 
		     :initform (error "Must initialize confidence-level."))
   (confidence-min :initarg :confidence-min :accessor confidence-min 
		   :initform (error "Must initialize confidence-min."))
   (confidence-max :initarg :confidence-max :accessor confidence-max 
		   :initform (error "Must initialize confidence-max."))
   (relative-error :initarg :relative-error :accessor relative-error 
		   :initform (error "Must initialize relative-error."))
   (absolute-error :initarg :absolute-error :accessor absolute-error 
		   :initform (error "Must initialize absolute-error."))
   (max-counts :initarg :max-counts :accessor max-counts 
	       :initform (error "Must initialize max-counts."))
   (binned-data :initarg :binned-data :accessor binned-data 
		:initform (error "Must initialize binned-data."))))



(defun make-parameter-distribution (result p no-bins start end confidence-level &optional no-iterations)
  (declare (ignore no-iterations))
  (let+ (((&values binned-data median min max max-counts)
	  (bin-parameter-values result p
				:no-bins no-bins
				:start start :end end
				:confidence-level confidence-level))
	 (abs-error (/ (- max min) 2d0)))
    (make-instance 'parameter-distribution
		   :name p
		   :median median
		   :confidence-level confidence-level
		   :confidence-min min
		   :confidence-max max
		   :relative-error (/ abs-error median)
		   :absolute-error abs-error
		   :max-counts max-counts
		   :binned-data binned-data)))

(defmethod get-parameter-info ((optimized-params optimized-parameters) parameter-name)
  (let+ (((&slots parameter-infos) optimized-params))
    (alexandria:if-let (p (find parameter-name parameter-infos :key #'name))
      p (error 'ba:unknown-parameter :parameter-name parameter-name))))







