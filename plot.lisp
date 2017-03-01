(in-package #:bayesian-analysis)


;;; api

(defgeneric plot-result-model (result &key))
(defgeneric plot-iteration-values (result &key (params-to-plot) (start) end (every)))
(defgeneric plot-data (data &key))


;;; implementation

(defun %get-mgl-data-for-data-1d-w/errors (data)
  (let+ (((&slots independent-parameters dependent-parameters
		  error-parameters) data)
	 (xs (slot-value data (first independent-parameters)))
	 (ys (slot-value data (first dependent-parameters)))
	 (ss (slot-value data (first error-parameters))))
    (iter
      (for x in-sequence xs)
      (for y in-sequence ys)
      (for s in-sequence ss)
      (collect (list x y s) into pd)
      (minimize x into min-x)
      (maximize x into max-x)
      (finally (return (values pd min-x max-x))))))



(defmethod plot-data ((data data) &key (style "lc 7 pt 7 lw 1.5 pw 1.5")
				       title xlabel ylabel)
  (unless (= 1 (no-independent-parameters data)
	     (no-dependent-parameters data) (no-error-parameters data))
    (error "plotting of given data not supported."))
  (labels ((cmd (fmt-str &rest args)
	       (apply #'format t fmt-str args)
	     (mgl-gnuplot:command (apply #'format nil fmt-str args)))
	   (iff (obj else) (if obj obj else)))
    (let+ (((x) (iff xlabel (independent-parameters data)))
	   ((y) (iff ylabel (dependent-parameters data)))
	   (title (iff title (format nil "~a_i[~a_i]" y x)))
	   ((&values plot-data) (%get-mgl-data-for-data-1d-w/errors data))
	   (options (format nil "using 1:2:3 with errorbars title '~a' ~a" title style)))
      (cmd "set xlabel '~a'" x)
      (cmd "set ylabel '~a'" y)
      (mgl-gnuplot:plot*
       (list
	(mgl-gnuplot:data* plot-data options))))))


(defmethod plot-result-model ((result solved-parameters) &key (no-steps 1000))
  (let+ (((&slots model data algorithm-result) result)
	 ((&slots input-model) algorithm-result)
	 (fun (model-function input-model)))
    (unless (= 1 (no-independent-parameters data)
	       (no-dependent-parameters data) (no-error-parameters data))
      (error "plotting of given data not supported."))
    (let+ (((&values plot-data min-x max-x) (%get-mgl-data-for-data-1d-w/errors data))
	   ((&values model-input-data model-results-data)
	    (iter
	      (for x from min-x to max-x by (/ (- max-x min-x) no-steps))
	      (collect (list x (funcall fun x input-model)) into id)
	      (collect (list x (funcall fun x model)) into rd)
	      (finally (return (values id rd))))))
      (mgl-gnuplot:plot*
       (list
	(mgl-gnuplot:data* plot-data "with errorbars title ''")
	(mgl-gnuplot:data* model-input-data "with lines lt 3 lw 0.3 lc 0 title 'input model'")
	(mgl-gnuplot:data* model-results-data "with lines lw 1.5 lc 7 title 'result model'"))))))


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








