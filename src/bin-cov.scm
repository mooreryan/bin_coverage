(declare (uses srfi-1
               srfi-13
               posix
               utils))
(declare (uses bc-utils))

(define VERSION "0.1.2")

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
  (system* "~a depth -aa ~a > ~a" SAMTOOLS bam-fname depth-fname))

(define (read-depth-file fname names-to-keep)

  (define (keep-node? node names-to-keep)
    (hash-table-exists? names-to-keep node))

  (define (parse-depth-line line)
    (let* ((line-split (string-split line "\t"))
           (node (first line-split))
           (pos (string->number (second line-split)))
           (cov (third line-split))
           (new-contig (= 1 pos)))
      (values node pos cov new-contig)))

  (define (keep-new-contig? new-contig node names-to-keep)
    (and new-contig (keep-node? node names-to-keep)))

  (define output-err-msg
    (string-join (list "No positions to graph."
                       "Do the names in the contig names file"
                       "match those in the bam/sam file?")
                 " "))

  (define (validate-and-return graph-lines start-posns names)
    (abort-if-empty graph-lines output-err-msg)
    (abort-if-empty start-posns output-err-msg)
    (abort-if-empty names output-err-msg)
    (values graph-lines start-posns names))

  (let ((port (open-input-file fname)))
    (define (rdf-iter graph-lines idx start-posns names keep)
      (let ((line (read-line port)))
        (if (eof-object? line)
            (validate-and-return graph-lines start-posns names)
            (let-values (((node pos cov new-contig) (parse-depth-line line)))
              (cond ((keep-new-contig? new-contig node names-to-keep)
                     (set! keep #t))
                    (new-contig (set! keep #f))
                    (else 'pass))
              (cond ((and new-contig keep)
                     (rdf-iter (cons (list idx cov) graph-lines)
                               (add1 idx)
                               (cons idx start-posns)
                               (cons node names)
                               keep))
                    (keep
                     (rdf-iter (cons (list idx cov) graph-lines)
                               (add1 idx)
                               start-posns
                               names
                               keep))
                    (else
                     (rdf-iter graph-lines
                               idx
                               start-posns
                               names
                               keep)))))))
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
  (system* "Rscript '~a'" rscript-fname))


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

(if (file-exists? depth-fname)
    (fprintf (current-error-port)
             "INFO -- Depth file (~a) exists, will use it~%"
             depth-fname)
    (begin
      (fprintf (current-error-port)
               "LOG -- Running samtools depth~%")
      (samtools-depth bam-fname depth-fname)))

(fprintf (current-error-port)
         "LOG -- Reading depth file~%")
(define-values
  (graph-lines start-posns names)
  (read-depth-file depth-fname names-to-keep))

(print-log "Printing graph lines")
(print-graph-lines graph-lines graph-lines-fname)

(print-log "Printing start posns lines")
(print-start-posns start-posns start-posns-fname)

(print-log "Printing contig names")
(print-start-posns names names-fname)

(print-log "Writing R script")
(write-rscript rscript-fname)

(print-log "Running rscript")
(run-rscript rscript-fname)

;; (print-log "Cleaning up")
;; (delete-file depth-fname)
;; (delete-file graph-lines-fname)
;; (delete-file start-posns-fname)
;; (delete-file names-fname)
;; (delete-file rscript-fname)
