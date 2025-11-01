
! This module was automatically created with Mathematica!
MODULE Hamiltonian
  IMPLICIT NONE
  !
  COMPLEX(8), PARAMETER :: ci = CMPLX(0,1,KIND=8)
  REAL(8), PARAMETER :: pi = 3.14159265358979323846264_8
  ! REAL(8), PARAMETER :: E = 2.71828182845904523536028_8
  REAL(8), PARAMETER :: s2 = Sqrt(2._8), s3 = Sqrt(3._8)
  !
  CONTAINS
  ! 
  SUBROUTINE hamilton(H,kx,ky,kz,P)
  COMPLEX(8), DIMENSION(1,1) :: H
  REAL(8), DIMENSION(4) :: P
  REAL(8) :: kx, ky, kz
  !
  REAL(8) :: ckx, cky, ckz
  REAL(8) :: ckxd2, ckyd2, ckzd2
  REAL(8) :: skx, sky, skz
  REAL(8) :: skxd2, skyd2, skzd2
  !
  ckx = COS(kx*Pi)
  cky = COS(ky*Pi)
  ckz = COS(kz*Pi)
  ! 
  ckxd2 = COS((kx*Pi)/2._8)
  ckyd2 = COS((ky*Pi)/2._8)
  ckzd2 = COS((kz*Pi)/2._8)
  !
  skx = SIN(kx*Pi)
  sky = SIN(ky*Pi)
  skz = SIN(kz*Pi)
  ! 
  skxd2 = SIN((kx*Pi)/2._8)
  skyd2 = SIN((ky*Pi)/2._8)
  skzd2 = SIN((kz*Pi)/2._8)
  !
  H(1,1) = P(1) + 4._8*(ckx*cky + (ckx + cky)*ckz)*P(2) + &
           2._8*(Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi) + &   
           Cos(2._8*kz*Pi))*P(3) +8._8*ckx*cky*ckz*P(4)	
  END SUBROUTINE
END MODULE
