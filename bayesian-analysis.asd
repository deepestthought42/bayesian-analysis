;;;; bayesian-analysis.asd

(asdf:defsystem #:bayesian-analysis
  :description "Package to specify a model and calculate posterior distributions as well
  as odds ratios using Bayesian statistics."
  :author "Renee Klawitter <klawitterrenee@gmail.com>"
  :license "Apache 2.0"
  :version "0.0.1"
  :depends-on (#:alexandria
               #:iterate
               #:let-plus
	       #:gsl-cffi
	       #:mgl-gnuplot
	       #:math-utils
	       #:cffi-nlopt
	       #:lla
	       #:verbose)
  :serial t
  :components ((:file "package")
	       (:file "bayesian-analysis")
	       (:file "utils")
	       (:file "init")
	       (:file "algorithm")
	       (:file "data")
	       (:file "priors")
	       (:file "sampling")
	       (:file "likelihood")
	       (:file "model-function")
               (:file "model")
	       (:file "metropolis-hastings")
	       (:file "levenberg-marquardt")
	       (:file "integration")
	       (:file "nlopt")
	       (:file "fisher-information")
	       (:file "odds")
	       (:file "plot")))

