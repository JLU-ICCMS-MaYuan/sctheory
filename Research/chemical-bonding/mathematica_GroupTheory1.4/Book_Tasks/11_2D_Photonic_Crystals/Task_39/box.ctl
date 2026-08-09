; Compute band structure for a square lattice of dielectric rods
; in air.

; Define various parameters with define-param so that they are
; settable from the command-line (with mpb <param>=<value>):
(define-param r 0.2) ; radius of the rods
(define-param eps 4.0) ; dielectric constant
(define-param k-interp 50) ; number of k points to interpolate
;(define-param phi 0.0) ; rotation angle for object

(define GaAs (make dielectric (epsilon eps)))

(set! geometry-lattice (make lattice (size 1 1 no-size))) ; 2d cell

(set! geometry 
      (list
       (make block
         (material GaAs)
         (center 0.5 0.5) (size .5 .5 no-size)
         (e1 1 0 0) (e2 0 1 0) (e3 0 0 1)
       )
      )
)
      
(define Gamma (vector3 0 0 0))
(define X (vector3 0.5 0 0))
(define M (vector3 0.5 0.5 0))
(set! k-points (interpolate k-interp (list Gamma X M Gamma)))

(set-param! resolution 32)
(set-param! num-bands 20)

; Compute the TE and TM bands.  Wrap in the (begin-time message ...)
; construct from libctl so that we report the total elapsed time:
(begin-time
 "total time for both TE and TM bands: "
 (run-te)
 (run-tm))

(display-eigensolver-stats)
