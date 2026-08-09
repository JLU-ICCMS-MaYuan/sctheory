!****m* src/tb_data.f90
!
! NAME
!  TB_Data.f90
!
! AUTHOR
!  W. Hergert
!
! USE 
!  The module is used to get the ab initio band structure data 
!  into the subroutines of fitting procedures
! 
! SOURCE
!
MODULE tbdata
!
USE typedef
IMPLICIT NONE
INTEGER(I4B)                             :: nkp    ! number of k-points
INTEGER(I4B)                             :: nbands ! maximum number of bands
INTEGER(I4B)                             :: nparm  ! number of TB-parameters
INTEGER(I4B)                             :: nenerg ! number of energies for fit
INTEGER(I4B)                             :: hdim   ! dimension Hamiltonian
INTEGER(I4B)                             :: nbu    ! lower band index from Hamiltonian used in fit
INTEGER(I4B)                             :: nbo    ! upper band index
!
REAL(DP), AllOCATABLE, DIMENSION(:,:)    :: kmesh  ! k-points used for fit
REAL(DP), AllOCATABLE, DIMENSION(:,:)    :: eabin  ! ab initio bandstructure 
REAL(DP), AllOCATABLE, DIMENSION(:,:)    :: eweight! weight for bands (VASP)
REAL(DP), AllOCATABLE, DIMENSION(:)      :: kweight! weight for k-points (VASP)
ReAL(DP), ALLOCATABLE, DIMENSION(:)      :: weight ! weights for fit
INTEGER(I4B), AllOCATABLE, DIMENSION(:)  :: sel    ! upper band number
INTEGER(I4B), AllOCATABLE, DIMENSION(:)  :: seu    ! lower band number
!
REAL(DP), AllOCATABLE, DIMENSION(:)      :: efit   ! energies for fit
!
! Data for GA
!
REAL(DP)                                 :: xmin,xmax ! Rescaling
INTEGER                                  :: seed      ! start random numbers
LOGICAL                                  :: lscal     ! type rescaling
REAL(DP), ALLOCATABLE, DIMENSION(:)      :: scal      ! Rescaling
END MODULE tbdata
!***
