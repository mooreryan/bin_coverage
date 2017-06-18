(declare (unit bc-utils))

(define (abort-if-empty lst . msgs)
  (if (zero? (length lst))
      (begin
        (fprintf (current-error-port)
                 "ERROR -- ~a~%"
                 (string-join msgs " "))
        (exit 1))))

(define (check-file fname)
  (if (not (file-exists? fname))
      (begin (fprintf (current-error-port)
                      "ERROR -- File '~a' does not exist"
                      fname)
             (exit 1))))

(define (print-log msg)
  (fprintf (current-error-port)
           "LOG -- ~a~%"
           msg))

(define (try-mkdir path)
  (if (directory-exists? outdir)
      (begin
        (fprintf (current-error-port)
                 "ERROR -- the directory '~a' already exists~%"
                 outdir)
        (exit 1))
      (create-directory outdir)))
