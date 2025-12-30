;
; set details for the calculation
;
(set-param! resolution 32)
(set-param! num-bands 8)
;
; compute the TE and TM bands. time measurement used 
;

(run-tm)

(display-eigensolver-stats)
