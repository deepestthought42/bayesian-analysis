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
	       #:verbose)
  :serial t
  :components ((:file "package")
	       (:file "init")
	       (:file "algorithm")
	       (:file "data")
	       (:file "priors")
	       (:file "sampling")
	       (:file "likelihood")
	       (:file "model-function")
               (:file "model")
	       (:file "metropolis-hastings")
	       (:file "odds")
	       (:file "bayesian-analysis")
	       (:file "plot")))

