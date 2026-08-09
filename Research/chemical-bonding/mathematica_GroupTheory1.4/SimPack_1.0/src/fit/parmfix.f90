!****m* src/tparmfix.f90
!
! NAME
!  parmfix.f90
!
! AUTHOR
!  W. Hergert
!
! USE
!  Select from a parameter set only those which are not fixed for fitting
!  Return complete data set after fitting.
!
! INPUT
!  x    - restricted set for fit
!  ps   - full set including fixed parameters
!  isel - = -1 restrict set
!         = +1 complete set
!
! COMMON /pfix/
!  ifp  - parameter numbers to be fixed
!  vfix - values of the fixed parameters
!  nump - number of parameters
!  nfix - number of fixed parameters
! 
! SOURCE
!
  SUBROUTINE parmfix(x,ps,isel)
  USE typedef
  REAL(DP)         :: vfix, x(1),ps(1)
  INTEGER(I4B)     :: ifp,nump,nfix,isel,l,k,i
  LOGICAL          :: lfix
  COMMON /pfix/ vfix(40),ifp(40),nump,nfix,lfix
! WRITE(*,*) nump,nfix,lfix
! WRITE(*,*) vfix
! WRITE(*,*) ifp
!
  IF(isel.LT.0) THEN
!
! select restricted set for fitting
! 
    l=1
    k=1
    DO i=1,nump
       IF(i.EQ.ABS(ifp(k)))THEN
         k=k+1
       ELSE
         x(l)=ps(i)        
         l=l+1
       ENDIF
    ENDDO    
  ELSEIF(isel.GT.0) THEN
!
! construct complete set
!    
    l=1
    k=1
    DO i=1,nump
       IF(i.EQ.ifp(k)) THEN      
          ps(i)=vfix(k)
          k=k+1
       ELSE
         ps(i)=x(l)
         l=l+1
       ENDIF  
    ENDDO
   
  ENDIF
   RETURN  
!
END SUBROUTINE parmfix
