!****m* src/tb_simplex.f90
!
! NAME
!  tb_simplex.f90
!
! AUTHOR
!  W. Hergert
!
! USE
!  SUbroutine for DNLQ1
!  Calculates the deviation of the band structure 
!  model from the ab initio data
!
! INPUT
!  N   - number of parameters in the model
!  M   - number of equations, i.e. ab initio energies
!  X   - paramter vector
!  FX  - vector containing the deviations
!  IER - error indicator (=0 all is o.k.)
!
! SOURCE
!
     DOUBLE PRECISION FUNCTION tb_simplex(x)
     USE typedef
     USE tbdata
     USE Hamiltonian
!
! COMMON Block to fix parameters of the set
!
     COMMON /pfix/ vfix(40),ifp(40),nump,nfix,lfix
     REAL(DP)      :: vfix
     INTEGER(I4B)  :: ifp,nump,nfix
     LOGICAL       :: lfix
!         
!    IMPLICIT NONE
!     INTEGER, PARAMETER        :: rk=kind(1.0D+00)
!
     REAL(DP), INTENT(IN)      :: x(nparm)
!     
     REAL(DP)                                  :: kx,ky,kz,diff,fx    
     INTEGER                                   :: kp,k,l,info,lwork
     CHARACTER(LEN=1)                          :: jobz,uplo
     REAL(DP), ALLOCATABLE, DIMENSION(:)       :: ev,rwork,work,xp
     COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:)  :: hamop
!
     uplo="U"
     jobz="N"
     lwork=3*hdim
     ALLOCATE(ev(hdim))
     ALLOCATE(work(3*hdim))
     ALLOCATE(rwork(3*hdim-2))
     ALLOCATE(hamop(hdim,hdim))  
!    
! Full set of parameters  
!
     ALLOCATE(xp(nump))
     IF(lfix) THEN
        CALL parmfix(x,xp,+1)
     ELSE
        xp=x
     ENDIF  
!     
!
! Loop over k-points
! 
     k=0
     fx=0.D0
! 
     DO kp=1,nkp
        kx=kmesh(1,kp)
        ky=kmesh(2,kp)
        kz=kmesh(3,kp)
        CALL hamilton(hamop,kx,ky,kz,xp)
!
        CALL zheev(jobz,uplo,hdim,hamop,hdim,ev,work,lwork,rwork,info)
!   
        IF(info.NE.0) THEN
          WRITE(ip,*) 'zheev   ',kp,info
         STOP
        ENDIF
        DO l=1,sel(kp)
           k=k+1            
           diff=(efit(k)-ev(l))*weight(l)
           fx=fx+diff*diff
        ENDDO
     ENDDO
     fx=SQRT(fx)/nkp/hdim
     tb_simplex=fx
     RETURN
!
     DEALLOCATE(ev,work,rwork,hamop,xp)
     END FUNCTION tb_simplex
!***	
