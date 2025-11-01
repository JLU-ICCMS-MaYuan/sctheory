
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
  COMPLEX(8), DIMENSION(9,9) :: H
  REAL(8), DIMENSION(23) :: P
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
  H(1,1) = 4._8*(ckx*cky + (ckx + cky)*ckz)*P(1) + 2._8*(Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi) + Cos(2._8*kz*Pi))*P(11) + P(21)
  H(2,1) = -(ci*(2._8*s2*(ckx + ckz)*P(2)*sky + 2._8*P(12)*Sin(2._8*ky*Pi)))
  H(3,1) = -(ci*(2._8*s2*(ckx + cky)*P(2)*skz + 2._8*P(12)*Sin(2._8*kz*Pi)))
  H(4,1) = -(ci*(2._8*s2*(cky + ckz)*P(2)*skx + 2._8*P(12)*Sin(2._8*kx*Pi)))
  H(5,1) = -2._8*s3*P(5)*skx*sky
  H(6,1) = -2._8*s3*P(5)*sky*skz
  H(7,1) = (-2._8*ckx*cky + (ckx + cky)*ckz)*P(5) + (-Cos(2._8*kx*Pi) - Cos(2._8*ky*Pi) + 2._8*Cos(2._8*kz*Pi))*P(15)
  H(8,1) = -2._8*s3*P(5)*skx*skz
  H(9,1) = s3*(ckx - cky)*ckz*P(5) + s3*(Cos(2._8*kx*Pi) - Cos(2._8*ky*Pi))*P(15)
  H(1,2) = 2._8*s2*ci*(ckx + ckz)*P(2)*sky + 2._8*ci*P(12)*Sin(2._8*ky*Pi)
  H(2,2) = 2._8*(cky*(ckx + ckz)*P(3) + (ckx*cky + (2._8*ckx + cky)*ckz)*P(4)) + 2._8*(Cos(2._8*ky*Pi)*P(13) + &
           (Cos(2._8*kx*Pi) + Cos(2._8*kz*Pi))*P(14)) + P(22)
  H(3,2) = -2._8*P(3)*sky*skz + 2._8*P(4)*sky*skz
  H(4,2) = -2._8*P(3)*skx*sky + 2._8*P(4)*skx*sky
  H(5,2) = -(ci*(s2*(s3*cky*P(6) + 2._8*ckz*P(7))*skx + 2._8*P(17)*Sin(2._8*kx*Pi)))
  H(6,2) = -(ci*(s2*(s3*cky*P(6) + 2._8*ckx*P(7))*skz + 2._8*P(17)*Sin(2._8*kz*Pi)))
  H(7,2) = -(ci*(-(((2._8*ckx*P(6) + ckz*(-P(6) + 2._8*s3*P(7)))*sky)/s2) - P(16)*Sin(2._8*ky*Pi)))
  H(8,2) = 0._8
  H(9,2) = -(ci*(-(((s3*ckz*P(6) + 2._8*(2*ckx + ckz)*P(7))*sky)/s2) - s3*P(16)*Sin(2._8*ky*Pi)))
  H(1,3) = 2._8*s2*ci*(ckx + cky)*P(2)*skz + 2._8*ci*P(12)*Sin(2._8*kz*Pi)
  H(2,3) = 2._8*(-P(3) + P(4))*sky*skz
  H(3,3) = 4._8*ckx*cky*P(4) + 2._8*(ckx + cky)*ckz*(P(3) + P(4)) + 2._8*(Cos(2._8*kz*Pi)*P(13) + &
           (Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi))*P(14)) + P(22)
  H(4,3) = -2._8*P(3)*skx*skz + 2._8*P(4)*skx*skz
  H(5,3) = 0._8
  H(6,3) = -(ci*(s2*(s3*ckz*P(6) + 2._8*ckx*P(7))*sky + 2._8*P(17)*Sin(2._8*ky*Pi)))
  H(7,3) = -(ci*(((ckx + cky)*(P(6) + 2._8*s3*P(7))*skz)/s2 + 2._8*P(16)*Sin(2._8*kz*Pi)))
  H(8,3) = -(ci*(s2*(s3*ckz*P(6) + 2._8*cky*P(7))*skx + 2._8*P(17)*Sin(2._8*kx*Pi)))
  H(9,3) = ci*(-(Sqrt(1.5_8)*ckx*P(6)*skz) + Sqrt(1.5_8)*cky*P(6)*skz + s2*ckx*P(7)*skz - s2*cky*P(7)*skz)
  H(1,4) = 2._8*s2*ci*(cky + ckz)*P(2)*skx + 2._8*ci*P(12)*Sin(2._8*kx*Pi)
  H(2,4) = 2._8*(-P(3) + P(4))*skx*sky
  H(3,4) = 2._8*(-P(3) + P(4))*skx*skz
  H(4,4) = 4._8*cky*ckz*P(4) + 2._8*ckx*(cky + ckz)*(P(3) + P(4)) + 2._8*(Cos(2._8*kx*Pi)*P(13) + &
          (Cos(2._8*ky*Pi) + Cos(2._8*kz*Pi))*P(14)) + P(22)
  H(5,4) = -(ci*(s2*(s3*ckx*P(6) + 2._8*ckz*P(7))*sky + 2._8*P(17)*Sin(2._8*ky*Pi)))
  H(6,4) = 0._8
  H(7,4) = -(ci*(-(((2._8*cky*P(6) + ckz*(-P(6) + 2._8*s3*P(7)))*skx)/s2) - P(16)*Sin(2._8*kx*Pi)))
  H(8,4) = -(ci*(s2*(s3*ckx*P(6) + 2._8*cky*P(7))*skz + 2._8*P(17)*Sin(2._8*kz*Pi)))
  H(9,4) = -(ci*(((s3*ckz*P(6) + 2._8*(2*cky + ckz)*P(7))*skx)/s2 + s3*P(16)*Sin(2._8*kx*Pi)))
  H(1,5) = -2._8*s3*P(5)*skx*sky
  H(2,5) = s2*ci*(s3*cky*P(6) + 2._8*ckz*P(7))*skx + 2._8*ci*P(17)*Sin(2._8*kx*Pi)
  H(3,5) = 0._8
  H(4,5) = s2*ci*(s3*ckx*P(6) + 2._8*ckz*P(7))*sky + 2._8*ci*P(17)*Sin(2._8*ky*Pi)
  H(5,5) = ckx*cky*(3._8*P(8) + P(10)) + 2._8*(ckx + cky)*ckz*(P(9) + P(10)) + 2._8*((Cos(2._8*kx*Pi) + & 
           Cos(2._8*ky*Pi))*P(19) + Cos(2._8*kz*Pi)*P(20)) + P(23)
  H(6,5) = -2._8*P(9)*skx*skz + 2._8*P(10)*skx*skz
  H(7,5) = s3*P(8)*skx*sky - s3*P(10)*skx*sky
  H(8,5) = -2._8*P(9)*sky*skz + 2._8*P(10)*sky*skz
  H(9,5) = 0._8
  H(1,6) = -2._8*s3*P(5)*sky*skz
  H(2,6) = s2*ci*(s3*cky*P(6) + 2._8*ckx*P(7))*skz + 2._8*ci*P(17)*Sin(2._8*kz*Pi)
  H(3,6) = s2*ci*(s3*ckz*P(6) + 2._8*ckx*P(7))*sky + 2._8*ci*P(17)*Sin(2._8*ky*Pi)
  H(4,6) = 0._8
  H(5,6) = 2._8*(-P(9) + P(10))*skx*skz
  H(6,6) = 2._8*ckx*cky*(P(9) + P(10)) + ckz*(cky*(3._8*P(8) + P(10)) + 2._8*ckx*(P(9) + P(10))) + &
           2._8*((Cos(2._8*ky*Pi) + Cos(2._8*kz*Pi))*P(19) + Cos(2._8*kx*Pi)*P(20)) + P(23)
  H(7,6) = -(s3*P(8)*sky*skz)/2._8 + (s3*P(10)*sky*skz)/2._8
  H(8,6) = -2._8*P(9)*skx*sky + 2._8*P(10)*skx*sky
  H(9,6) = (3._8*P(8)*sky*skz)/2._8 - (3._8*P(10)*sky*skz)/2._8
  H(1,7) = (-2._8*ckx*cky + (ckx + cky)*ckz)*P(5) + (-Cos(2._8*kx*Pi) - Cos(2._8*ky*Pi) + 2._8*Cos(2._8*kz*Pi))*P(15)
  H(2,7) = -((ci*(2._8*ckx*P(6) + ckz*(-P(6) + 2._8*s3*P(7)))*sky)/s2) - ci*P(16)*Sin(2._8*ky*Pi)
  H(3,7) = (ci*(ckx + cky)*(P(6) + 2._8*s3*P(7))*skz)/s2 + 2._8*ci*P(16)*Sin(2._8*kz*Pi)
  H(4,7) = -((ci*(2._8*cky*P(6) + ckz*(-P(6) + 2._8*s3*P(7)))*skx)/s2) - ci*P(16)*Sin(2._8*kx*Pi)
  H(5,7) = s3*(P(8) - P(10))*skx*sky
  H(6,7) = (s3*(-P(8) + P(10))*sky*skz)/2._8
  H(7,7) = (cky*ckz*(P(8) + 12._8*P(9) + 3._8*P(10)) + ckx*(4._8*cky*(P(8) + 3._8*P(10)) + ckz*(P(8) + 12._8*P(9) + &
           3._8*P(10))))/4._8 + 2._8*Cos(2._8*kz*Pi)*P(18) + ((Cos(2._8*kx*Pi) + &
           Cos(2._8*ky*Pi))*(P(18) + 3._8*P(20)))/2._8 + P(23)
  H(8,7) = -(s3*P(8)*skx*skz)/2._8 + (s3*P(10)*skx*skz)/2._8
  H(9,7) = -(s3*(-ckx + cky)*ckz*(P(8) - 4._8*P(9) + 3._8*P(10)))/4._8 - (s3*(-Cos(2._8*kx*Pi) +  &
           Cos(2._8*ky*Pi))*(-P(18) + P(20)))/2._8
  H(1,8) = -2._8*s3*P(5)*skx*skz
  H(2,8) = 0._8
  H(3,8) = s2*ci*(s3*ckz*P(6) + 2._8*cky*P(7))*skx + 2._8*ci*P(17)*Sin(2._8*kx*Pi)
  H(4,8) = s2*ci*(s3*ckx*P(6) + 2._8*cky*P(7))*skz + 2._8*ci*P(17)*Sin(2._8*kz*Pi)
  H(5,8) = 2._8*(-P(9) + P(10))*sky*skz
  H(6,8) = 2._8*(-P(9) + P(10))*skx*sky
  H(7,8) = (s3*(-P(8) + P(10))*skx*skz)/2._8
  H(8,8) = 2._8*ckx*cky*(P(9) + P(10)) + ckz*(ckx*(3._8*P(8) + P(10)) + 2._8*cky*(P(9) + P(10))) + &
           2._8*((Cos(2._8*kx*Pi) + Cos(2._8*kz*Pi))*P(19) + Cos(2._8*ky*Pi)*P(20)) + P(23)
  H(9,8) = (-3._8*P(8)*skx*skz)/2._8 + (3._8*P(10)*skx*skz)/2._8
  H(1,9) = s3*(ckx - cky)*ckz*P(5) + s3*(Cos(2._8*kx*Pi) - Cos(2._8*ky*Pi))*P(15)
  H(2,9) = -((ci*(s3*ckz*P(6) + 2._8*(2*ckx + ckz)*P(7))*sky)/s2) - s3*ci*P(16)*Sin(2._8*ky*Pi)
  H(3,9) = (ci*(-ckx + cky)*(-(s3*P(6)) + 2._8*P(7))*skz)/s2
  H(4,9) = (ci*(s3*ckz*P(6) + 2._8*(2*cky + ckz)*P(7))*skx)/s2 + s3*ci*P(16)*Sin(2._8*kx*Pi)
  H(5,9) = 0._8
  H(6,9) = (3._8*(P(8) - P(10))*sky*skz)/2._8
  H(7,9) = -(s3*(-ckx + cky)*ckz*(P(8) - 4._8*P(9) + 3._8*P(10)))/4._8 - (s3*(-Cos(2._8*kx*Pi) + &
           Cos(2._8*ky*Pi))*(-P(18) + P(20)))/2._8
  H(8,9) = (3._8*(-P(8) + P(10))*skx*skz)/2._8
  H(9,9) = 4._8*ckx*cky*P(9) + ((ckx + cky)*ckz*(3._8*P(8) + 4._8*P(9)))/4._8 + (9._8*(ckx + cky)*ckz*P(10))/4._8 + &
           (3._8*(Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi))*P(18) + (Cos(2._8*kx*Pi) + Cos(2._8*ky*Pi) + &
           4._8*Cos(2._8*kz*Pi))*P(20))/2._8 + P(23)
  END SUBROUTINE
END MODULE