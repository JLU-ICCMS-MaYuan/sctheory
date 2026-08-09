!****m* src/dqsf.f90
!
! NAME
!  dqsf.f90
!
! AUTHOR
!  W. Hergert
!
! USE 
!  Integration by means of Simpsons formula
! 
! INPUTS
!
!  H    - step width of aequidistant data points
!  Y    - curve to integrate
!  Z    - integrated curve 
!  NDIM - number of data points (>3) 
!
! SOURCE
!
      SUBROUTINE DQSF(H,Y,Z,NDIM)
!
USE typedef
! 
REAL(DP), INTENT(IN)    :: Y(NDIM),H
INTEGER,  INTENT(IN)    :: NDIM
REAL(DP), INTENT(OUT)   :: Z(NDIM)

INTEGER                 :: i,j,m
!
      IF(NDIM.LE.3) THEN
         WRITE(il,1) NDIM
1        FORMAT(1X,'*** ERROR SUBROUTINE DQSF *** NDIM=',I2)
         STOP
      ENDIF   
      Z(1)=Y(1)
      DO I=2,NDIM
         IF(I.EQ.2.AND.NDIM.EQ.3) Z(I)=H*(5.D0*Y(1)/4.D0+2.D0*Y(2)-Y(3)/4.D0)/3.D0
         IF(I.EQ.2.AND.NDIM.GT.3) Z(I)=3.D0*H/8.D0*(Y(1)+3.D0*Y(2)+3.D0*Y(3)+  &
                                       Y(4))-H/3.D0*(Y(2)+4.D0*Y(3)+Y(4))
         IF(I.EQ.4) Z(I)=3.D0*H/8.D0*(Y(1)+3.D0*Y(2)+3.D0*Y(3)+Y(4))
         IF(I.EQ.2.OR.I.EQ.4) GOTO 200
         M=I-2
         IF(MOD(I,2).EQ.0) GOTO 300
         Z(I)=0.D0
         DO J=1,M,2
            Z(I)=Z(I)+H/3.D0*(Y(J)+4.D0*Y(J+1)+Y(J+2))
         ENDDO
         GOTO 200
300      Z(I)=H/24.D0*(8.*(Y(1)+Y(6))+31.D0*(Y(2)+Y(5))+21.D0*(Y(3)+Y(4)))
         IF(I.EQ.6) GOTO 200
         DO J=6,M,2
            Z(I)=Z(I)+H/3.D0*(Y(J)+4.D0*Y(J+1)+Y(J+2))
         ENDDO
200   CONTINUE
      ENDDO
      RETURN
      END
