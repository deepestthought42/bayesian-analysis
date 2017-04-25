(in-package #:bayesian-analysis)


(defparameter *debug-function* nil
  "When non-nil *DEBUG-FUNCTION* is assumed to be a function with
  lamba-list: (CATEGORY FORMAT-CONTROL &REST FORMAT-ARGS) to do debug
  output.")

(defparameter *lev/marq-error-absolute* 1d-6)
(defparameter *lev/marq-error-relative* 1d-6)
(defparameter *lev/marq-no-max-iterations 100)


(defun default-debug-out (level cat fmt-ctrl &rest fmt-args)
  (apply #'v:log level cat fmt-ctrl fmt-args))

(defun debug-out (level category fmt &rest params)
  (apply *debug-function* level category fmt params))








