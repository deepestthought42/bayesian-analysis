;;;; bayesian-analysis.asd

(asdf:defsystem #:bayesian-analysis
  :description "Describe bayesian-analysis here"
  :author "Your Name <your.name@example.com>"
  :license "Specify license here"
  :depends-on (#:alexandria
               #:iterate
               #:let-plus
	       #:gsl-cffi
	       #:mgl-gnuplot
	       #:math-utils
	       #:cffi-nlopt
	       #:verbose)
  :serial t
  :components ((:file "package")
	       (:file "bayesian-analysis")
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
	       (:file "fisher-information")
	       (:file "odds")
	       (:file "plot")))

