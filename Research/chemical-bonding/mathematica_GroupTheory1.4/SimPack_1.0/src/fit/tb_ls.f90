!****m* src/tb_deviate.f90
!
! NAME
!  tb_deviate.f90
!
! AUTHOR
!  W. Hergert
!
! USE
!  SUbroutine for DNLQ1
!  Calculates the deviation of the band structure 
!  model from the ab initio data
!
!  The module EWP is used if MKL is not available.
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
     SUBROUTINE tb_ls(N,X,M,FX,F1,IER)
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
!
     REAL(DP), INTENT(IN)      :: x(n)
     REAL(DP), INTENT(OUT)     :: fx(m)
     INTEGER(I4B), INTENT(IN)  :: n,m
     INTEGER(I4B), INTENT(OUT) :: ier
     EXTERNAL F1
!     
     REAL(DP)                                    :: kx,ky,kz    
     INTEGER                                     :: kp,k,l,info,lwork
     CHARACTER(LEN=1)                            :: jobz,uplo
     REAL(DP), ALLOCATABLE, DIMENSION(:)         :: ev,rwork,work,xp
     COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:)    :: hamop
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
! Loop over k-points
! 
     k=0
     DO kp=1,nkp
        kx=kmesh(1,kp)
        ky=kmesh(2,kp)
        kz=kmesh(3,kp)
        CALL hamilton(hamop,kx,ky,kz,xp)
! 
        CALL zheev(jobz,uplo,hdim,hamop,hdim,ev,work,lwork,rwork,info)
        IF(info.NE.0) THEN
           WRITE(ip,*) kp,info
           STOP
        ENDIF
!     
        DO l=1,sel(kp)
           k=k+1
           fx(k)=(efit(k)-ev(l))*weight(l)
        ENDDO
     ENDDO
!
     IER=0
     DEALLOCATE(ev,work,rwork,hamop,xp)
     END SUBROUTINE tb_ls
!
     SUBROUTINE F1(N)  
     USE typedef  
     IMPLICIT NONE
     INTEGER(I4B) :: N
     END SUBROUTINE F1 
!***
