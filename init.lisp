(in-package #:bayesian-analysis)


(defparameter *debug-function* nil)


(defun debug-out (category fmt &rest params)
  (declare (ignore category))
  (apply *debug-function* fmt params))




