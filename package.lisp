;;;; package.lisp

(defpackage #:bayesian-analysis
  (:use #:cl #:let-plus #:iterate)
  (:export
   #:approximate
   #:model
   #:unknown-prior-type
   #:unknown-sampling-type
   #:mcmc-iterations
   #:initialize-accumulator
   #:initialize-likelihood
   #:data
   #:data-with-sigma
   #:initialize-from-source
   #:y
   #:no-data-points
   #:sigma)
  (:nicknames ba))

