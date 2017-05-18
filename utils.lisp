(in-package #:bayesian-analysis)


(defstruct val max min)


(defun get-min/max (sequence &key (max most-negative-double-float)
				  (min most-positive-double-float))
  (let ((ret (reduce #'(lambda (val no)
			 (setf (val-max val) (max no (val-max val))
			       (val-min val) (min no (val-min val)))
			 val)
		     sequence
		     :initial-value (make-val :max max :min min))))
    (values (val-min ret) (val-max ret))))



(defun sumlogexp (sequence)
  "Returns ln[∑exp[xᵢ]], where the xᵢ are member of the sequence given
in SEQUENCE. It is assuming the elementes of SEQUENCE to be of type
DOUBLE-FLOAT and minimizes the losses due to floating point under/overflow."
  (let+ (((&values min max) (get-min/max sequence))
	 (n (length sequence))
	 ;; dividing by n here as all entries could be bigger/smaller than
	 ;; than what most/least positve would allow
	 (max-log (log (/ most-positive-double-float n)))
	 (min-log (log (* least-positive-double-float n)))
	 (a (cond
	      ;; no shift needed
	      ((and (< max max-log)
		    (> min min-log))
	       0d0)
	      ;; as infinity is worse than zero -- since sbcl at least
	      ;; return 0d0 in that case -- always shift to what
	      ;; should never be infinite. This also covers the case
	      ;; where (and (> max max-log) (< min min-log)).
	      ((> max max-log)
	       (- max max-log))
	      ;; if shifting by (- min min-log) keeps the maximum
	      ;; below max-log, do it -- otherwise, shift so max =
	      ;; max-log
	      ((< min min-log)
	       (if (<= (+ max (- min-log min)) max-log)
		   (- min-log min)
		   (- max-log max)))
	      (t (error 'program-error)))))
    (log (+ a (reduce #'(lambda (x1 x2) (+ x2 (exp (- x1 a)))) sequence)))))



















