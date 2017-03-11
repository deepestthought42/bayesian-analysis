(in-package #:bayesian-analysis)


;;; api

(defgeneric plot-result-model (parameter-result &key))
(defgeneric plot-iteration-values (mcmc-result &key (params-to-plot) (start) end (every)))
(defgeneric plot-data (data &key))
(defgeneric plot-likelihood (mcmc-result &key (start) end (every)))

;;; implementation


(defun %get-mgl-data (data x-slot other-slots)
  (let+ (((&values min-x max-x)
	  (iter
	    (for x in-sequence (slot-value data x-slot))
	    (minimize x into min-x)
	    (maximize x into max-x)
	    (finally (return (values min-x max-x)))))
	 (len (length (slot-value data x-slot)))
	 (offset (+ min-x (/ (- max-x min-x) 2))))
    (values
     (iter
       (for i from 0 below len)
       (collect
	   (append (list (- (aref (slot-value data x-slot) i) offset))
		   (iter
		     (for slot in other-slots)
		     (collect (aref (slot-value data slot) i))))))
     min-x max-x offset)))

(defun %get-gnuplot-labels-for-1x/1y-and-title (x-slot y-slot title xlabel ylabel offset)
  (labels ((cmd (fmt-str &rest args)
	       (apply #'format t fmt-str args)
	     (mgl-gnuplot:command (apply #'format nil fmt-str args)))
	   (iff (obj else) (if obj obj else)))
    (let+ ((x (iff xlabel x-slot))
	   (xx (format nil "~a - ~f" x offset))
	   (y (iff ylabel y-slot)))
      (values xx y
	      (iff title (format nil "~a_i[~a_i - ~f]" y x offset))))))


(defun get-plot-mgl-plot-data (data independent-parameter other-parameters
			       plot-type style title xlabel ylabel)
  (labels ((cmd (fmt-str &rest args)
	     (apply #'format t fmt-str args)
	     (mgl-gnuplot:command (apply #'format nil fmt-str args))))
    (let+ (((&values plot-data min-x max-x offset)
	    (%get-mgl-data data
			   independent-parameter
			   other-parameters))
	   ((&values x y title) (%get-gnuplot-labels-for-1x/1y-and-title
				 independent-parameter (first other-parameters) title xlabel ylabel offset))
	   (options (format nil "with ~a  ~a title '~a'"
			    plot-type style title)))
      (values
       (mgl-gnuplot:data* plot-data options) x y min-x max-x offset))))



(defun get-plot-mgl-data-depending-on-type (data style title xlabel ylabel)
  (labels ((f (slot) (first (slot-value data slot))))
    (cond
      ((= 1 (no-independent-parameters data)
	  (no-dependent-parameters data)
	  (no-error-parameters data))
       (get-plot-mgl-plot-data data
			       (f 'independent-parameters)
			       (list (f 'dependent-parameters)
				     (f 'error-parameters))
			       "errorbars" style
			       title xlabel ylabel))
      ((and (= 1 (no-independent-parameters data) (no-dependent-parameters data))
	    (= 0 (no-error-parameters data)))
       (get-plot-mgl-plot-data data
			       (f 'independent-parameters)
			       (list (f 'dependent-parameters))
			       "points" style title xlabel ylabel))
      ((and (= 2 (no-independent-parameters data))
	    (= (no-error-parameters data)
	       (no-dependent-parameters data)
	       0))
;       (mgl-gnuplot:command "set style fill transparent solid 0.35 noborder")
       (get-plot-mgl-plot-data data
			       (f 'independent-parameters)
			       (cdr (slot-value data 'independent-parameters))
			       "points solid " style title xlabel ylabel))
      (t (error "plotting of given data not supported.")))))

(defmethod plot-data ((data data) &key (style "lc 7 pt 7 lw 1 pw 1")
				       title xlabel ylabel)
  (labels ((cmd (fmt-str &rest args)
	     (apply #'format t fmt-str args)
	     (mgl-gnuplot:command (apply #'format nil fmt-str args))))
    (let+ (((&values plot-data x y min-x max-x offset)
	    (get-plot-mgl-data-depending-on-type data style title xlabel ylabel)))
      (unless (= min-x max-x)
	(cmd "set xrange [~f:~f]" (- min-x offset) (- max-x offset)))
      (cmd "set xlabel '~a'" x)
      (cmd "set ylabel '~a'" y)
      (mgl-gnuplot:plot* (list plot-data)))))





(defmethod plot-result-model ((result solved-parameters)
			      &key (no-steps 1000)
				   (style-options/data "pt 7")
				   (style-options/input "with lines lt 3 lw 0.3 lc 0 title 'input input-model'")
				   (style-options/result "with lines lw 1.5 lc 7 title 'result input-model'"))
  (let+ (((&slots model data algorithm-result) result)
	 ((&slots input-model) algorithm-result)
	 (fun (model-function input-model)))
    (let+ (((&values plot-data x-label y-label min-x max-x offset)
	    (get-plot-mgl-data-depending-on-type data style-options/data "" nil nil))
	   ((&values model-input-data model-results-data)
	    (if (= min-x max-x)
		(values `((,(- min-x offset) ,(funcall fun min-x input-model)))
			`((,(- min-x offset) ,(funcall fun min-x model))))
		(iter
		  (for x from min-x to max-x by (/ (- max-x min-x) no-steps))
		  (collect (list (- x offset) (funcall fun x input-model)) into id)
		  (collect (list (- x offset) (funcall fun x model)) into rd)
		  (finally (return (values id rd)))))))
      (labels ((cmd (fmt-str &rest args)
		 (apply #'format t fmt-str args)
		 (mgl-gnuplot:command (apply #'format nil fmt-str args))))
	(cmd "set xrange [~f:~f]" (- min-x offset) (- max-x offset))
	(cmd "set xlabel '~a'" x-label )
	(cmd "set ylabel '~a'" y-label))
      (mgl-gnuplot:plot*
       (list
	plot-data
	(mgl-gnuplot:data* model-input-data style-options/input)
	(mgl-gnuplot:data* model-results-data style-options/result))))))




(defmethod plot-iteration-values ((result mcmc-parameter-result)
				  &key (params-to-plot) (start 0) (end) (every 1))
  (let+ (((&slots iteration-accumulator) result)
	 ((&slots marginalized-parameters parameter-array no-iterations) iteration-accumulator)
	 (end (if end end no-iterations)))
    (alexandria:if-let (diff (set-difference params-to-plot marginalized-parameters))
      (error "Unknown parameters: ~{~a~^, ~}" diff))
    (mgl-gnuplot:plot*
     (iter 
       (for p in (if params-to-plot params-to-plot marginalized-parameters))
       (for index-param = (position p marginalized-parameters))
       (collect
	   (mgl-gnuplot:data*
	    (iter
	      (for i from start below end)
	      (if (= (mod i every) 0)
	       (collect (list i (aref parameter-array index-param i)))))
	    (format nil "with steps title '~a'" p)))))))

(defmethod plot-parameter-distribution ((result solved-parameters) parameter)
  (let+ (((&slots parameter-infos) result)
	 (param (find parameter parameter-infos :key #'name))
	 ((&slots binned-data median confidence-min confidence-max max-counts) param))
    (labels ((cmd (fmt-str &rest args)
	       (apply #'format t fmt-str args)
	       (mgl-gnuplot:command (apply #'format nil fmt-str args))))
      (iter
	(for d in binned-data)
	(incf (car d) (- median)))
      (cmd "set style fill solid 0.1 noborder")
      (cmd "set arrow from ~,10f,0 to ~,10f,~,10f nohead front lt 1 lw 2 lc 7" 0 0 (* 1.01 max-counts))
      (mgl-gnuplot:plot*
       (list
	(mgl-gnuplot:data* (iter
			     (for (x c) in binned-data)
			     (if (<= (- confidence-min median) x (- confidence-max median))
				 (collect (list x c))))
			   (format nil "u 1:2 with boxes lc 7 title ''"))
	(mgl-gnuplot:data* binned-data (format nil "u 1:2 with histeps lc 0 title '~a'" parameter)))))))





(defmethod plot-likelihood ((result mcmc-parameter-result) &key (start 0) (end) (every 1))
  (let+ (((&slots iteration-accumulator input-model data) result)
	 (new-model (copy-object input-model))
	 (likelihood (initialize-likelihood new-model data))
	 ((&slots marginalized-parameters parameter-array no-iterations) iteration-accumulator)
	 (end (if end end no-iterations)))
    (mgl-gnuplot:plot*
     (list
      (mgl-gnuplot:data*
       (iter outer
	 (for i from start below end)
	 (when (= (mod i every) 0)
	   (collect
	       (iter
		 (for p in marginalized-parameters)
		 (for index-param = (position p marginalized-parameters))
		 (setf (slot-value new-model p) (aref parameter-array index-param i))
		 (finally (return (list i (funcall (varying/log-of-likelihood likelihood)))))))))
       (format nil "with steps title 'likelihood'"))))))


