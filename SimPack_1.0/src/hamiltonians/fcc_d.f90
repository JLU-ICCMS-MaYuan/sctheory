
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
  COMPLEX(8), DIMENSION(5,5) :: H
  REAL(8), DIMENSION(7) :: P
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
  H(1,1) = ckx*cky*(3._8*P(1) + P(3)) + 2._8*(ckx + cky)*ckz*(P(2) + P(3)) + 2._8*((Cos(2._8*kx*Pi) &
           + Cos(2._8*ky*Pi))*P(5) + Cos(2._8*kz*Pi)*P(6)) + P(7)
  H(2,1) = -2._8*P(2)*skx*skz + 2._8*P(3)*skx*skz
  H(3,1) = s3*P(1)*skx*sky - s3*P(3)*skx*sky
  H(4,1) = -2._8*P(2)*sky*skz + 2._8*P(3)*sky*skz
  H(5,1) = 0._8
  H(1,2) = 2._8*(-P(2) + P(3))*skx*skz
  H(2,2) = cky*ckz*(3._8*P(1) + P(3)) + 2._8*ckx*(cky + ckz)*(P(2) + P(3)) + 2._8*((Cos(2._8*ky*Pi) &
           + Cos(2._8*kz*Pi))*P(5) + Cos(2._8*kx*Pi)*P(6)) + P(7)
  H(3,2) = -(s3*P(1)*sky*skz)/2._8 + (s3*P(3)*sky*skz)/2._8
  H(4,2) = -2._8*P(2)*skx*sky + 2._8*P(3)*skx*sky
  H(5,2) = (3._8*P(1)*sky*skz)/2._8 - (3._8*P(3)*sky*skz)/2._8
  H(1,3) = s3*(P(1) - P(3))*skx*sky
  H(2,3) = (s3*(-P(1) + P(3))*sky*skz)/2._8
  H(3,3) = (cky*ckz*(P(1) + 12._8*P(2) + 3._8*P(3)) + ckx*(4._8*cky*(P(1) + 3._8*P(3)) + ckz*(P(1) &
            +12._8*P(2) + 3._8*P(3))))/4._8 + (4._8*Cos(2._8*kz*Pi)*P(4) + (Cos(2._8*kx*Pi) &
            + Cos(2._8*ky*Pi))*(P(4) + 3._8*P(6)))/2._8 + P(7)
  H(4,3) = -(s3*P(1)*skx*skz)/2._8 + (s3*P(3)*skx*skz)/2._8
  H(5,3) = -(s3*(-ckx + cky)*ckz*(P(1) - 4._8*P(2) + 3._8*P(3)))/4._8 - (s3*(-Cos(2._8*kx*Pi) &
           + Cos(2._8*ky*Pi))*(-P(4) + P(6)))/2._8
  H(1,4) = 2._8*(-P(2) + P(3))*sky*skz
  H(2,4) = 2._8*(-P(2) + P(3))*skx*sky
  H(3,4) = (s3*(-P(1) + P(3))*skx*skz)/2._8
  H(4,4) = 2._8*ckx*cky*(P(2) + P(3)) + ckz*(ckx*(3._8*P(1) + P(3)) + 2._8*cky*(P(2) + P(3)))  &
           +2._8*((Cos(2._8*kx*Pi) + Cos(2._8*kz*Pi))*P(5) + Cos(2._8*ky*Pi)*P(6)) + P(7)
  H(5,4) = (-3._8*P(1)*skx*skz)/2._8 + (3._8*P(3)*skx*skz)/2._8
  H(1,5) = 0._8
  H(2,5) = (3._8*(P(1) - P(3))*sky*skz)/2._8
  H(3,5) = -(s3*(-ckx + cky)*ckz*(P(1) - 4._8*P(2) + 3._8*P(3)))/4._8 - (s3*(-Cos(2._8*kx*Pi) &
           + Cos(2._8*ky*Pi))*(-P(4) + P(6)))/2._8
  H(4,5) = (3._8*(-P(1) + P(3))*skx*skz)/2._8
  H(5,5) = (16._8*ckx*cky*P(2) + (ckx + cky)*ckz*(3._8*P(1) + 4._8*P(2)) &
           + 9._8*(ckx + cky)*ckz*P(3))/4._8 + (3._8*(Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi))*P(4) &
           + (Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi) + 4._8*Cos(2._8*kz*Pi))*P(6))/2._8 + P(7)
  END SUBROUTINE
END MODULE