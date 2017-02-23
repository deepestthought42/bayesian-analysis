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



(defgeneric bin-parameter-values (result parameter bin-size &key start end))


(defun %calculate-confidance (array no-iterations confidance-level)
  (let+ ((sorted (sort (copy-seq array) #'> :key #'second)))
    (iter
      (for (x counts) in-sequence sorted)
      (sum counts into summed)
      (maximize x into max)
      (minimize x into min)
      (if (>= (/ summed no-iterations) confidance-level)
	  (return (values min max))))))





(defmethod bin-parameter-values ((result mcmc-parameter-result) parameter bin-width
				 &key (start 0) end (confidence-level 0.69))
  (let+ (((&slots iteration-accumulator) result)
	 (bin-width (coerce bin-width 'double-float))
	 ((&slots marginalized-parameters parameter-array
		  no-iterations) iteration-accumulator)
	 (end (if (not end) no-iterations end))
	 (pos (position parameter marginalized-parameters)))
    (if (not pos)
	(error "Parameter: ~a was not marginalized over." parameter))
    (let+ (((&values vals min max)
	    (iter
	      (for i from start below end)
	      (for val = (* bin-width (floor (aref parameter-array pos i) bin-width)))
	      (maximize val into max)
	      (minimize val into min)
	      (collect val into vals)
	      (finally (return (values vals min max)))))
	   (no-bins (round (/ (- max min) bin-width)))
	   (binned (make-array (1+ no-bins)
				  :initial-contents
				  (iter
				    (for i from 0 to no-bins)
				    (for x = (+ min (* (+ i 0.5) bin-width)))
				    (collect (list x 0))))))
      ;; binning
      (iter
	(for v in vals)
	(incf (cadr (aref binned (round (- v min) bin-width)))))
      ;; processing
      (iter
	(with median = (/ (- max min) 2))
	(with first-time = t)
	(with median-index = 0)
	(with counts = 0)
	(for (x c) in-sequence binned with-index i)
	(maximize c into max-counts)
	(incf counts c)
	(if (and first-time (>= counts (/ no-iterations 2)))
	    (setf first-time nil
		  median (car (aref binned i))
		  median-index i))
	(finally
	 (let+ (((&values min max)
		 (%calculate-confidance binned no-iterations confidence-level)))
	   (return (values (map 'list #'identity binned)
			   median min max max-counts))))))))










