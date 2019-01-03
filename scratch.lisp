(in-package #:bayesian-analysis)







(define-data-class 1d-data (x "x") y err
    (object (source t))
  (setf x (make-array 4 :initial-contents '(-1d0 0d0 1d0 2d0)
			:element-type 'double-float)
	y (make-array 4 :initial-contents '(2.99d0 3.1d0 3.1d0 2.99d0)
			:element-type 'double-float)
	err (make-array 4 :initial-contents '(0.1d0 0.1d0 0.2d0 0.1d0)
			  :element-type 'double-float)))


;(initialize-from-source '1d-data t)



(define-bayesian-model (quadratic 1d-data)
    ((a :default 0.5 :min -1 :max 1 :prior-type :uniform :sample-sigma 0.1d0)
     (b :prior-type :uniform :default -0.5)
     (c :prior-type :uniform :default 2 :min 2 :max 4 :sample-sigma 0.1d0))
    (:d_i=f_i+gaussian_error_i_unequal_sigma)
    ((x) (+ (* a x) (* b x x) c)))

(define-bayesian-model (linear 1d-data)
    ((a :default 1 :min -1 :max 1 :prior-type :uniform :sample-sigma 0.1d0)
     (b :prior-type :uniform :default 2 :min 2 :max 4 :sample-sigma 0.1d0))
    (:d_i=f_i+gaussian_error_i)
    ((x)
      (+ (* a x) b)))


(defparameter *test-model* (make-instance 'quadratic))


(labels ((cmd (fmt-str &rest args)
	   (mgl-gnuplot:command (apply #'format nil fmt-str args))))
  (mgl-gnuplot:with-session ()
    (cmd "reset")
    (cmd "set terminal x11 enhanced font 'Georgia,8' dashed")
    (plot-iteration-values
     (optimize (make-instance 'metropolis-hastings :no-iterations 50000)
			   (make-instance 'quadratic)
			   (initialize-from-source '1d-data t))
     :every 1 :end 5000) 
    (cmd "unset output")))


(labels ((cmd (fmt-str &rest args)
	   (mgl-gnuplot:command (apply #'format nil fmt-str args))))
  (mgl-gnuplot:with-session ()
    (cmd "reset")
q    (cmd "set terminal x11 enhanced font 'Georgia,12' dashed")
    (plot-parameter-distribution
     (get-parameter-results
      (optimize (make-instance 'metropolis-hastings :no-iterations 500000)
			    (make-instance 'quadratic :b-bin-width 0.01 :a-bin-width 0.01 :c-bin-width 0.004)
			    (initialize-from-source '1d-data t))
      :confidence-level 0.92
      :start 200)
     'c ) 
    (cmd "unset output")))



(labels ((cmd (fmt-str &rest args)
	   (mgl-gnuplot:command (apply #'format nil fmt-str args))))
  (mgl-gnuplot:with-session ()
    (cmd "reset")
    (cmd "set terminal wxt enhanced font 'Georgia,12' dashed")
    (plot-result
     (get-parameter-results
      (optimize (make-instance 'metropolis-hastings :no-iterations 100000)
		(make-instance 'quadratic :b-bin-width 0.001 :a-bin-width 0.001)
		(initialize-from-source '1d-data t))
      :start 200)) 
    (cmd "unset output")))


(labels ((cmd (fmt-str &rest args)
	   (mgl-gnuplot:command (apply #'format nil fmt-str args))))
  (mgl-gnuplot:with-session ()
    (cmd "reset")
    (cmd "set terminal wxt enhanced font 'Georgia,8' dashed")
    (plot-data (initialize-from-source '1d-data t))
    (cmd "unset output")))




(let+ ((r1 (get-parameter-results
	    (optimize
	     (make-instance 'metropolis-hastings :no-iterations 100000)
	     (make-instance 'linear :b-bin-width 0.005 :a-bin-width 0.005)
	     (initialize-from-source '1d-data t)) :start 200))
       (r2 (get-parameter-results
	    (optimize
	     (make-instance 'metropolis-hastings :no-iterations 100000)
	     (make-instance 'quadratic :b-bin-width 0.005 :a-bin-width 0.005)
	     (initialize-from-source '1d-data t)) :start 200))
       (m1 (model r1))
       (m2 (model r2)))
;  (incf (a m1) 0.01)
  (labels ((cmd (fmt-str &rest args)
	     (mgl-gnuplot:command (apply #'format nil fmt-str args))))
    (mgl-gnuplot:with-session ()
      (cmd "reset")
      (cmd "set terminal wxt enhanced font 'Georgia,12'")
      (plot-result r2) 
      (cmd "unset output")))
  (calculate-odds-ratio-1/2 m1 m2 (initialize-from-source '1d-data t)))
