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
     FUNCTION  tb_ga(np,x) RESULT(fs)
     USE typedef
     USE tbdata
     USE Hamiltonian
!     
     IMPLICIT NONE
!
     INTEGER, INTENT(IN)                    :: np
     REAL,    INTENT(IN)                    :: x(np)
     REAL(DP)                               :: fx
     REAL                                   :: fs

!     
     REAL(DP)                                :: kx,ky,kz,diff    
     INTEGER                                 :: kp,k,l,info,lwork,i
     CHARACTER(LEN=1)                        :: jobz,uplo
     REAL(DP), ALLOCATABLE, DIMENSION(:)     :: ev,rwork,work,xp
     COMPLEX(DP), ALLOCATABLE,DIMENSION(:,:) :: hamop  
!
! Rescaling of parameters for calculation 
!
     ALLOCATE(xp(np))
        DO i=1,np
           IF(lscal) THEN
               xp(i)=x(i)/scal(i)
           ELSE    
              xp(i)=x(i)*(xmax-xmin)+xmin
           ENDIF             
        ENDDO 
!    WRITE(ip,*) 'drin'
!    WRITE(ip,*) iscal
!    WRITE(ip,*)  xp 
!       
     uplo="U"
     jobz="N"
     lwork=3*hdim
     ALLOCATE(ev(hdim))
     ALLOCATE(work(3*hdim))
     ALLOCATE(rwork(3*hdim-2))
     ALLOCATE(hamop(hdim,hdim))  
  
!
! Loop over k-points
! 
     k=0
     fx=0.0_DP
     DO kp=1,nkp
!
!  a factor of 2 is temporarily included, due to a mistake in the exported
!  Hamiltonian    
!
        kx=kmesh(1,kp) !*2.0_DP
        ky=kmesh(2,kp) !*2.0_DP
        kz=kmesh(3,kp) !*2.0_DP
 !       WRITE(ip,*) kx,ky,kz
        CALL hamilton(hamop,kx,ky,kz,xp)
        CALL zheev(jobz,uplo,hdim,hamop,hdim,ev,work,lwork,rwork,info)
!        WRITE(ip,*) 'ev' ,ev
        IF(info.NE.0) THEN
           WRITE(ip,*) kp,info
           STOP
        ENDIF
        DO l=1,sel(kp)
           k=k+1            
!           write(ip,*) efit(k),'  ',ev(l)
           diff=efit(k)-ev(l)
           fx=fx+diff*diff
        ENDDO
     ENDDO
!
! calculate fitness
!     
     fs=1./SNGL(fx)
!
     DEALLOCATE(ev,work,rwork,hamop,xp)
!
     RETURN     
     END FUNCTION tb_ga
!***
