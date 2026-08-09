!****m* src/tb_abinitio.f90
!
! NAME
!  tb_abinitio.f90
!
! AUTHOR
!  W. Hergert
!
! USE 
!   The subroutine is used to read data from ab initio codes for 
!   tight-binding parameter fitting. ab initio data are stored in 
!   module tbdata.
!
! INPUT
!   file - file name of the data
!   code - selects the source of data 
!       STANDARD - a standard format mainly used for testing
!           VASP - VASP calculations
!  
! SOURCE
!
SUBROUTINE abinitio(file,code,verbose)
!
USE typedef
USE tbdata
IMPLICIT NONE
!
INTEGER(I4B)       :: verbose,i,k,numb,l,nb
REAL(DP)           :: ew,fac 
CHARACTER(LEN= 8)  :: code
CHARACTER(LEN=20)  :: file 
CHARACTER(LEN=80)  :: dummy
!
OPEN(8,FILE=file,FORM='formatted',STATUS='old')
!
!-----------------------------------------------------------------------
!     Selected method - STANDARD
!-----------------------------------------------------------------------    
!
IF(ADJUSTR(code).EQ.'STANDARD') THEN
  IF(verbose.GE.2) THEN
     WRITE(ip,*) '*** start tb_abinitio'
     WRITE(ip,*)
     WRITE(ip,*) '  Read data in Standard format'
  ENDIF
  READ(8,*) nkp
  READ(8,*) nbands
!
  ALLOCATE(kmesh(3,nkp))
  ALLOCATE(eabin(nkp,nbands))
  ALLOCATE(sel(nkp))
!
  DO i=1,nkp
     READ(8,*) sel(i),kmesh(1,i),kmesh(2,i),kmesh(3,i),(eabin(i,k),k=1,nbands)    
  ENDDO  
!
! total number of energies for fit
!
  nenerg=SUM(sel)  
  ALLOCATE(efit(nenerg))
!
  IF(verbose.GE.2) THEN
    WRITE(ip,10) nkp
    WRITE(ip,20) nbands
    WRITE(ip,40) nenerg
    WRITE(ip,*)
    WRITE(ip,*) '  Energies for fitting'
    WRITE(ip,*) '  k-point              select  energies '
    WRITE(ip,*) '*************************************************************'
    DO i=1,nkp
       WRITE(ip,30) kmesh(1,i),kmesh(2,i),kmesh(3,i), sel(i),(eabin(i,k),k=1,nbands)
    ENDDO   
  ENDIF
!
! vector of energy values
!
  k=0
  DO i=1,nkp
     DO l=1,sel(i)
        k=k+1
        efit(k)=eabin(i,l)
     ENDDO
  ENDDO
!
!-----------------------------------------------------------------------
!     Selected method - VASP
!-----------------------------------------------------------------------    
!
ELSEIF(ADJUSTR(code).EQ.'    VASP') THEN
 IF(verbose.GE.2) THEN
  WRITE(ip,*) '*** start tb_abinitio'
     WRITE(ip,*)
    WRITE(ip,*) 'Read VASP data'
      WRITE(ip,*)
 ENDIF
 DO i=1,5
    READ(8,*) dummy
    IF(verbose.GE.2) THEN
      WRITE(ip,'(A60)') dummy
    ENDIF
 ENDDO
 WRITE(ip,*)
 READ(8,*) numb,nkp,nb  
!
! check data
! 
 IF(nb.LT.nbands) THEN
   WRITE(ip,*) 'Error: VASP data do not contain eneough bands.'
   STOP
   ELSE
 ENDIF
!
! Allocate array.
!
 ALLOCATE(kmesh(3,nkp))
 ALLOCATE(eabin(nkp,nb))
 ALLOCATE(eweight(nkp,nb))
 ALLOCATE(kweight(nkp))
 IF(verbose.GE.1) THEN
    WRITE(ip,20) nb
    WRITE(ip,10) nkp
 ENDIF   
 DO k=1,nkp
    READ(8,*) kmesh(1,k),kmesh(2,k),kmesh(3,k),kweight(k)
    DO i=1,nb
       READ(8,*) l, eabin(k,i),eweight(k,i)
    ENDDO 
 ENDDO   
  If(verbose.EQ.2) THEN
    WRITE(ip,*)
    WRITE(ip,*) '  Energies for fitting (VASP)'
    WRITE(ip,*) '  k-point                 energies '
    WRITE(ip,*) '*************************************************************'
    DO i=1,nkp
       WRITE(ip,50) kmesh(1,i),kmesh(2,i),kmesh(3,i),(eabin(i,k),k=1,nb)
    ENDDO   
  ENDIF
!
!-----------------------------------------------------------------------
!     Select method unknown
!-----------------------------------------------------------------------    
!
ELSE
  WRITE(ip,*) 'Method not implemented'
  STOP 
END IF  
  IF(verbose.GE.2) THEN
     WRITE(ip,*)
          WRITE(ip,*) '*** end tb_abinitio'
           WRITE(ip,*)
            WRITE(ip,*)
  ENDIF
!
!***********************************************************************
!
10 FORMAT(3x,'Number of k-points         : nkp    = ',I3)
20 FORMAT(3x,'Number of bands            : nbands = ',I3)
40 FORMAT(3x,'Number of energies for fit : nenerg = ',I3)
30 FORMAT(3x,3(F6.4,1x),I5,2x,20(F8.4,1x))
50 FORMAT(3x,3(F6.4,1x),2x,10(F8.4,1x)/(26x,10(F8.4,1x)))
!
END SUBROUTINE abinitio
!***
