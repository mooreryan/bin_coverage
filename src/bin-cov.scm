(declare (uses srfi-1
               posix))

(define VERSION "0.1.0")

;; The keys of the hash are the contig names for faster lookup
(define (read-contig-names fname)
  (let ((port (open-input-file fname))
        (ht (make-hash-table string=? string-hash)))
    (define (rcn-iter)
      (let ((line (read-line port)))
        (if (eof-object? line)
            ht
            (begin (hash-table-set! ht line 1)
                   (rcn-iter)))))
    (rcn-iter)))

(define (samtools-depth bam-fname depth-fname)
  (system (sprintf "~a depth -aa ~a > ~a" SAMTOOLS bam-fname depth-fname)))

;; (define (keep-node? node names-to-keep)
;;   (member node names-to-keep))

(define (keep-node? node names-to-keep)
  (hash-table-exists? names-to-keep node))

(define (read-depth-file fname names-to-keep)
  (let ((port (open-input-file fname)))
    (define (rdf-iter graph-lines idx start-posns names keep)
      (let ((line (read-line port)))
        (if (eof-object? line)
            (values (reverse graph-lines)
                    (reverse start-posns)
                    (reverse names))
            (let* ((line-split (string-split line "\t"))
                   (node (first line-split))
                   (pos (string->number (second line-split)))
                   (cov (third line-split))
                   (new-contig (= 1 pos)))
              ;; Check to see if keep should flip
              (cond ((and new-contig (keep-node? node names-to-keep))
                     (set! keep #t))
                    (new-contig (set! keep #f))
                    (else 'pass))
              (if keep
                  (rdf-iter (cons (list idx cov) graph-lines)
                            (add1 idx)
                            (if new-contig
                                (cons idx start-posns)
                                start-posns)
                            (if new-contig
                                (cons node names)
                                names)
                            keep)
                  (rdf-iter graph-lines
                            idx
                            start-posns
                            names
                            keep))))))
    (rdf-iter '() 1 '() '() #f)))

(define (print-graph-lines graph-lines graph-lines-fname)
  (let ((port (open-output-file graph-lines-fname)))
    (fprintf port "position coverage~%")
    (for-each (lambda (lst)
                (fprintf port
                         "~a ~a~%"
                         (first lst)
                         (second lst)))
              graph-lines)))

(define (print-start-posns start-posns start-posns-fname)
  (let ((port (open-output-file start-posns-fname)))
    (fprintf port "start.pos~%")
    (for-each (lambda (item)
                (fprintf port "~a~%" item))
              start-posns)))

(define (write-rscript fname)
  (define rscript
    (string-join
     (list
      (sprintf "dat <- read.table(\"~a\", sep=\" \", header=T);" graph-lines-fname)
      (sprintf "start.posns <- read.table(\"~a\", sep=\" \", header=T);" start-posns-fname)
      "width.mult <- nrow(dat) / 50000"
      (sprintf "png(\"~a\", units = \"in\", res=120, width = max(1 * width.mult, 11), height=5);" cov-plot-fname)
      "plot(dat, type = \"l\", xlab=\"Position\", ylab = \"Coverage\");"
      "lapply(start.posns, FUN=function(pos){abline(v=pos, col=\"blue\")});"
      (sprintf "names <- read.table(\"~a\", header=T);" names-fname)
      "for(idx in 1:nrow(names)){text(cex = 0.75, srt=90, col=\"red\", pos=4, x=start.posns$start.pos[idx], y = max(dat$coverage)/2, names$start.pos[idx])};"
      "dev.off();")
     "\n"))
  (let ((port (open-output-file fname)))
    (fprintf port "~a~%" rscript)))

(define (run-rscript rscript-fname)
  (system (sprintf "Rscript '~a'" rscript-fname)))

(define (try-mkdir path)
  (if (directory-exists? outdir)
      (begin
        (fprintf (current-error-port)
                 "ERROR -- the directory '~a' already exists~%"
                 outdir)
        (exit 1))
      (create-directory outdir)))

(define (check-file fname)
  (if (not (file-exists? fname))
      (begin (fprintf (current-error-port)
                      "ERROR -- File '~a' does not exist"
                      fname)
             (exit 1))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(if (not (= 3 (length (command-line-arguments))))
    (begin (fprintf (current-error-port)
                    "VERSION: v~a~%USAGE: ~a /path/to/samtools recruitment.bam contig-names.txt~%"
                    VERSION
                    (program-name))
           (exit 1)))

(for-each check-file (command-line-arguments))

(define SAMTOOLS (first (command-line-arguments)))
(define bam-fname (second (command-line-arguments)))
(define bin-contigs-fname (third (command-line-arguments)))

(define depth-fname (sprintf "~a.bin-cov.samtools-depth-aa" bam-fname))
(define graph-lines-fname (sprintf "~a.graph-lines" depth-fname))
(define start-posns-fname (sprintf "~a.start-posns" depth-fname))
(define names-fname (sprintf "~a.names" depth-fname))
(define cov-plot-fname (sprintf "~a.cov-plot.png" graph-lines-fname))
(define rscript-fname (sprintf "~a.rscript.r" depth-fname))

(fprintf (current-error-port)
         "LOG -- Reading contigs fname~%")
(define names-to-keep (read-contig-names bin-contigs-fname))

(fprintf (current-error-port)
         "LOG -- Running samtools depth~%")
(samtools-depth bam-fname depth-fname)

(fprintf (current-error-port)
         "LOG -- Reading depth file~%")
(define-values (graph-lines start-posns names) (read-depth-file depth-fname names-to-keep))


(fprintf (current-error-port)
         "LOG -- Printing graph lines~%")
(print-graph-lines graph-lines graph-lines-fname)

(fprintf (current-error-port)
         "LOG -- Printing start posns lines~%")
(print-start-posns start-posns start-posns-fname)

(fprintf (current-error-port)
         "LOG -- Printing contig names~%")
(print-start-posns names names-fname)

(fprintf (current-error-port)
         "LOG -- Writing R script~%")
(write-rscript rscript-fname)

(fprintf (current-error-port)
         "LOG -- Running rscript~%")
(run-rscript rscript-fname)

(fprintf (current-error-port)
         "LOG -- Cleaning up~%")
;; (delete-file depth-fname)
;; (delete-file graph-lines-fname)
;; (delete-file start-posns-fname)
;; (delete-file names-fname)
;; (delete-file rscript-fname)
