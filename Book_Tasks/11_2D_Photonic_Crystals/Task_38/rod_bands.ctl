;-----------------------------------------------------------------------
; Group Theory in Solid state physics - Problem solving with Mathematica
;-----------------------------------------------------------------------
;
; This input file corresponds to the example section of the MPB tutorial
;
; Band structure for a square lattice of dielectric rods
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
; define special k-points and the path in k-space
; 
(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))
(define M (vector3 0.5 0.5 0))

(define-param k-interp 50) ; number of k points to interpolate

(set! k-points (interpolate k-interp (list Gamma X M Gamma)))
;
; set details for the calculation
;
(set-param! resolution 32)
(set-param! num-bands 8)
;
; compute the TE and TM bands. time measurement used 
;
(begin-time
  "total time for both TE and TM bands: "
 (run-te)
 (run-tm)
)

(display-eigensolver-stats)
