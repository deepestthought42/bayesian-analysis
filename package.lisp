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
   #:sigma
   #:f-of-x
   #:get-all-data-slots
   #:get-independent-parameters
   #:get-dependent-parameters
   #:get-error-parameters
   #:update-model-with-results
   #:parameter-result
   #:plot-result-model
   #:plot-iteration-values
   #:plot-data
   #:calculate-odds-ratio
   #:incongruent-data
   #:wrong-data-type
   #:wrong-number-of-arguments
   #:get-parameter-description-text
   #:unknown-parameter)
  (:nicknames ba))

