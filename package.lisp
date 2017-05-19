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
   #:unknown-parameter
   #:define-bayesian-model
   #:define-data-class
   #:metropolis-hastings
   #:find-optimum
   #:get-parameter-results
   #:plot-parameter-distribution
   #:*debug-function*
   #:default-debug-out
   #:plot-likelihood
   #:new-sample?
   #:cache
   #:get-1d-plot-function
   #:varying/log-of-likelihood
   #:constant/log-of-likelihood
   #:get-prior-for-parameter
   #:gsl-multifit-fdfsolver-alloc
   #:no-parameters-to-marginalize
   #:nlopt
   #:nlopt-result
   #:copy-object
   #:gaussian-prior
   #:mu
   #:*default-sampler*
   #:gaussian-sampler)
  (:nicknames ba))

