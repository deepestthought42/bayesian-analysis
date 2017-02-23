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
(defgeneric get-parameter-results (result &key start end))

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



(defclass solved-parameters ()
  ((model :initarg :model :accessor model 
	  :initform (error "Must initialize model."))
   (parameter-infos :initarg :parameter-infos :accessor parameter-infos 
		    :initform (error "Must initialize parameter-infos."))
   (algorithm-result :initarg :algorithm-result :accessor algorithm-result 
		     :initform (error "Must initialize algorithm-result."))
   (data :initarg :data :accessor data 
	 :initform (error "Must initialize data."))))


(defclass parameter-result ()
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
   (max-counts :initarg :max-counts :accessor max-counts 
	       :initform (error "Must initialize max-counts."))
   (binned-data :initarg :binned-data :accessor binned-data 
		:initform (error "Must initialize binned-data."))
   (bin-width :initarg :bin-width :accessor bin-width 
	      :initform (error "Must initialize bin-width."))))


(defmethod get-parameter-results ((result mcmc-parameter-result) &key (start 0) end (confidence-level 0.69))
  (let+ (((&slots iteration-accumulator input-model no-iterations data) result)
	 ((&slots model-parameters-to-marginalize) input-model)
	 (model (copy-object input-model))
	 (end (if end end no-iterations))
	 (param-infos (iter
			(for p in model-parameters-to-marginalize)
			(for bin-width = (slot-value model (w/suffix-slot-category :bin-width p)))
			(let+ (((&values binned-data median min max max-counts)
				(bin-parameter-values result p bin-width
						      :start start :end end
						      :confidence-level confidence-level)))
			  (setf (slot-value model p) median)
			  (collect (make-instance 'parameter-result
						  :name p
						  :median median
						  :confidence-level confidence-level
						  :confidence-min min
						  :confidence-max max
						  :max-counts max-counts
						  :binned-data binned-data
						  :bin-width bin-width))))))
    (make-instance 'solved-parameters
		   :algorithm-result result
		   :parameter-infos param-infos
		   :data data
		   :model model)))





