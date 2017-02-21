(in-package #:bayesian-analysis)






(define-data-class 1d-data x y err
    ((source) 0)
    (object (source t))
  (setf x (make-array 4 :initial-contents '(-1d0 0d0 1d0 2d0))
	y (make-array 4 :initial-contents '(2.99d0 3.1d0 3.1d0 2.99d0))
	err (make-array 4 :initial-contents '(0.1d0 0.1d0 0.2d0 0.1d0))))





(define-bayesian-model (quadratic 1d-data)
    ((a :default 2 :min -10 :max 10 :prior-type :uniform :sample-sigma 1d0)
     (b :prior-type :uniform :default 2 :min -10 :max 10 :sample-sigma 1d0)
     (c :prior-type :uniform :default 2 :min -10 :max 10 :sample-sigma 1d0))
    (:d_i=f_i+gaussian_error_i)
    ((x) (+ (* a x) (* b x x) c)))

(define-bayesian-model (linear 1d-data)
    ((a :default 1 :min -10 :max 10 :prior-type :uniform :sample-sigma 0.1d0)
     (b :prior-type :uniform :default 2 :min -10 :max 10 :sample-sigma 0.1d0))
    (:d_i=f_i+gaussian_error_i)
    ((x)
      (+ (* a x) b)))


(defparameter *test-model* (make-instance 'quadratic))



(let ((*debug-function* nil
			;; #'(lambda (fmt &rest args)
			;;     (apply #'format t fmt args)
			;;     (format t "~%"))
			))
  (solve-for-parameters (make-instance 'metropolis-hastings :no-iterations 50000)
			(make-instance 'linear)
			(initialize-from-source '1d-data t)))
