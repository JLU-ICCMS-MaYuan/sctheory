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
!
!REAL(DP), AllOCATABLE, DIMENSION(:,:)    :: kmesh  ! k-points used for fit
REAL(DP), AllOCATABLE, DIMENSION(:,:)    :: eabin  ! ab initio bandstructure 
REAL(DP), AllOCATABLE, DIMENSION(:,:)    :: eweight! weight for bands (VASP)
REAL(DP), AllOCATABLE, DIMENSION(:)      :: kweight! weight for k-points (VASP)
INTEGER(I4B), AllOCATABLE, DIMENSION(:)  :: sel    ! number bands per k-point
!
REAL(DP), AllOCATABLE, DIMENSION(:)      :: efit   ! energies for fit
!
! Data for GA
!
REAL(DP)                                 :: xmin,xmax ! Rescaling
INTEGER(I4B)                             :: seed      ! Seed random number generator
END MODULE tbdata
!***
