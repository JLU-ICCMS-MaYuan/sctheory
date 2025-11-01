;-----------------------------------------------------------------------
; Group Theory in Solid state physics - Problem solving with Mathematica
;-----------------------------------------------------------------------
;
; This input file corresponds to the example section of the MPB tutorial
;
; Calculate band structure at k-points for DOS calculation
;
; Define the real structure and permittivity profile
; 
(define-param r 0.2)       ; radius of the rods
(define-param eps 8.9)     ; dielectric constant

(define Alumina (make dielectric (epsilon eps)))

(set! geometry-lattice (make lattice (size 1 1 no-size))) ; 2d cell

(set! geometry 
      (list
       (make cylinder 
	 (material Alumina) 
	 (center 0 0) (radius r) (height infinity)
       )
      )
)
;
; include the GTPack output here!
; 

