!****g* src/typedef
! 
!  NAME
!   typedef
!  AUTHOR
!   W. Hergert
!  TYPE
!   MODULE
!  PACKAGE
!   to all packages
!  USAGE
!   definition of the types used in all programs
!  MODIFICATION HISTORY
!   August 2013   : first derived data types with respect to TB
!   December 2014 : additional derived data types
!
! SOURCE
!
MODULE typedef
!
! Symbolic names for INTEGERS and REALS as well as COMPLEX
!
INTEGER, PARAMETER :: IB = SELECTED_INT_KIND(18)
INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: SP  = SELECTED_REAL_KIND(6,37)
INTEGER, PARAMETER :: DP  = KIND(1.0D+00)
INTEGER, PARAMETER :: QP  = SELECTED_REAL_KIND(28,2465)
!
INTEGER, PARAMETER :: LGT = 1         ! ?
!
! parameters  fixed for input and output
!
INTEGER(I2B), PARAMETER :: ir=5       ! standard input
INTEGER(I2B), PARAMETER :: ip=6       ! standard output
INTEGER(I2B), PARAMETER :: il=20      ! unit for verbose log file
INTEGER(I2B), PARAMETER :: ierr=0     ! unit error -> screen
!
! mathematical constants
!
REAL(DP),    PARAMETER      :: delta  =EPSILON(1.D0)
REAL(DP),    PARAMETER      :: pi     =3.141592653589793238462643_DP
REAL(DP),    PARAMETER      :: sqr2   =1.414213562373095048801689_DP
COMPLEX(DP), PARAMETER      :: iu=(0.0_DP,1.0_DP) 
!
END MODULE typedef
!***
