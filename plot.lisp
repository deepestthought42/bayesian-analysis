(in-package #:bayesian-analysis)


;;; api

(defgeneric plot-result (optimization-result &key))
(defgeneric plot-result-models (input-model result-model data &key))

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


(defparameter *data-transparency* 0.2)
(defparameter *2d-data-circle-size* 1)

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
       (mgl-gnuplot:command (format nil "set style fill transparent solid ~,2f noborder" *data-transparency*))
       (mgl-gnuplot:command (format nil "set style circle radius ~,2f" *2d-data-circle-size*))
       (get-plot-mgl-plot-data data
			       (f 'independent-parameters)
			       (cdr (slot-value data 'independent-parameters))
			       "circles" style title xlabel ylabel))
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



(defmethod plot-result ((result optimized-parameters)
			&key (no-steps 1000)
			     (style-options/data "pt 7")
			     (style-options/input "with lines lt 3 lw 0.3 lc 0 title 'input input-model'")
			     (style-options/result "with lines lw 1.5 lc 7 title 'result input-model'")
			     (enclose-in-plot t)
			     (overwrite-xlabel))
  (let+ (((&slots model data algorithm-result) result)
	 ((&slots input-model) algorithm-result))
    (plot-result-models input-model model data
			:no-steps no-steps
			:style-options/data style-options/data
			:style-options/input style-options/input
			:style-options/result style-options/result
			:enclose-in-plot enclose-in-plot
			:overwrite-xlabel overwrite-xlabel)))

(defmethod plot-result-models ((input-model model) (result-model model) data 
			       &key (no-steps 1000)
				    (style-options/data "pt 7")
				    (style-options/input "with lines lt 3 lw 0.3 lc 0 title 'input input-model'")
				    (style-options/result "with lines lw 1.5 lc 7 title 'result input-model'")
				    (include-input-model t)
				    (include-data t)
				    (enclose-in-plot t)
				    (overwrite-xlabel))
  (let+ ((input-fun (ba:get-1d-plot-function input-model))
	 (result-fun (ba:get-1d-plot-function result-model))
	 ((&values plot-data x-label y-label min-x max-x offset)
	  (get-plot-mgl-data-depending-on-type data style-options/data "" nil nil))
	 ((&values model-input-data model-results-data)
	  (if (= min-x max-x)
	      (values `((,(- min-x offset) ,(funcall input-fun min-x)))
		      `((,(- min-x offset) ,(funcall result-fun min-x))))
	      (iter
		(for x from min-x to max-x by (/ (- max-x min-x) no-steps))
		(collect (list (- x offset) (funcall input-fun x)) into id)
		(collect (list (- x offset) (funcall result-fun x)) into rd)
		(finally (return (values id rd)))))))
    (labels ((cmd (fmt-str &rest args)
	       (apply #'format t fmt-str args)
	       (mgl-gnuplot:command (apply #'format nil fmt-str args))))
      (cmd "set xrange [~f:~f]" (- min-x offset) (- max-x offset))
      (cmd "set xlabel '~a'" (if overwrite-xlabel overwrite-xlabel x-label))
      (cmd "set ylabel '~a'" y-label))
    (let ((d `(,@(if include-data (list plot-data))
	       ,@(if include-input-model (list (mgl-gnuplot:data* model-input-data style-options/input)))
	       ,(mgl-gnuplot:data* model-results-data style-options/result))))
      (if enclose-in-plot (mgl-gnuplot:plot* d) d))))




(defmethod plot-iteration-values ((result mcmc-optimization-result)
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

(defmethod plot-parameter-distribution ((result optimized-parameters) parameter
					&key title (offset-around-median t)
					     (fn-on-x #'identity))
  (let+ (((&slots parameter-infos) result)
	 (param (alexandria:if-let (p (find parameter parameter-infos :key #'name))
		  p (error "Couldn't find parameter: ~a" parameter)))
	 ((&slots binned-data median confidence-min confidence-max max-counts) param)
	 (median (funcall fn-on-x median))
	 (confidence-min (funcall fn-on-x confidence-min))
	 (confidence-max (funcall fn-on-x confidence-max))
	 (offset (if offset-around-median median 0d0))
	 (binned-data
	        (iter
		  (for d in binned-data)
		  (collect (list (- (funcall fn-on-x (first d)) offset)
				 (second d))))))
    (labels ((cmd (fmt-str &rest args)
	       (apply #'format t fmt-str args)
	       (mgl-gnuplot:command (apply #'format nil fmt-str args))))
      (cmd "unset arrow")
      (cmd "set arrow from ~,10f,0 to ~,10f,~,10f nohead front lt 1 lw 2 lc 7"
	   (if offset-around-median 0d0 median)
	   (if offset-around-median 0d0 median)
	   (* 1.05 max-counts))
      (mgl-gnuplot:plot*
       (list
	(mgl-gnuplot:data* (iter
			     (for (x c) in binned-data)
			     (if (<=  (- confidence-min offset) x (- confidence-max offset))
				 (collect (list x c))))
			   (format nil "u 1:2 with boxes lc 7 fs solid 0.3 noborder title ''"))
	(mgl-gnuplot:data* binned-data (format nil "u 1:2 with histeps lc 0 title '~a'"
					       (cond
						 ((stringp title) title)
						 ((functionp title) (funcall title median))
						 (t (format nil "~a, median: ~,3f" parameter median))))))))))





(defmethod plot-likelihood ((result mcmc-optimization-result) &key (start 0) (end) (every 1))
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


