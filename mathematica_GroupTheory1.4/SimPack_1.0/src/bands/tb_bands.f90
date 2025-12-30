!****m* src/tb_bands.f90
!
! NAME
!  tb_bands.f90
! 
! AUTHOR
!  W. Hergert
!
! HISTORY
!  * 30.03.2022 first version
!  * 04.04.2022 input and output now similar to the DOS proram
!  * 27.04.2022 k-mesh generation from LMTO implmented
!  * 10.05.2022 check form
!  * 08.09.2022 check verbosity level
!  * 16.09.2022 changes band calculation for DOS
!  * 29.11.2022 restict bands for DOS calculation to one tetrahedron method
!
! USAGE 
!  Calculation of energy bands and eventually weight factors for a 
!  given set of k-points
!
! DESCRIPTION
!  This program is part of SimPack, a collection of simple FORTRAN 
!  routines to be used in connection with the Mathematica group theory 
!  package GTPack.
!
!  Usually calculations of band structures can be done directly in 
!  Mathematica. In some cases (large matrix of the Hamilonian and/or 
!  large number of k-points) this faster version can be used.
!
!  It is also helpful to test the FORTRAN export of the Hamiltonians 
!  from GTPack.
!
!  The program calculates also band structure data for DOS calculations.
!  Two methods are implemented:
!  GAUSS  - Gaussian smearing
!  TETRA  - tetrahedron method Dresden group (old code for cubic lattices)
! 
!  The k-vetors for GAUSS are mainly exported from GTPack, but also the 
!  k-meshes from TETRA can be used for the Gaussian smearing 
!  method in the DOS calculation.
!
! INPUTS
!  A set of tokens is defined to structure the input (CHARACTER(LEN=10)
!   
! File names :
!
!   'F_PROTOCOL'    - file name protocol            (out)
!   'F_BNDSTRUC'    - file name band structure data (out)
!   '  F_MATELM'    - file name matrix elements     (out)
!   'F_ETB_PARM'    - file name ETB parameters      (in)
!   'F_KVECTORS'    - file name k-vectors           (in)
!
! Other control data :
!
!   'PROBLMTEXT'    - headline for the calculation
!   '  DOS_CALC'    - method of DOS calculation
!   ' KMESH_GEN'    - data for generation of k-mesh for DOS 
!   'HAM_DIMENS'    - dimension of Hamiltonian
!   '   VERBOSE'    - definition oflevel of verbosity
!                      = 0 only tokens
!                      = 1 k-vectors used
!                      = 2 k-vectors and eigenvalues
!
! SOURCE
!
!
PROGRAM tb_bands
USE typedef
USE tbdata
USE Hamiltonian
USE kpoints_3D
!
IMPLICIT NONE
!
COMPLEX(DP), ALLOCATABLE, DIMENSION(:,:) :: ham
REAL(DP), ALLOCATABLE, DIMENSION(:,:)    :: bands,kmesh
REAL(DP), ALLOCATABLE, DIMENSION(:)      :: w,work,rwork,we,weight

!
INTEGER                                  :: i,j,k,m,lwork,info,imat,nbd,np,npart(20), &
                                            iverb,nteil,ist
CHARACTER(LEN=5)                         :: TLAT(15)
CHARACTER(LEN=6)                         :: dos
CHARACTER(LEN=1)                         :: jobz,uplo
CHARACTER(LEN=25)                        :: nam(50)
CHARACTER(LEN=60)                        :: problem,f_prot,f_band,f_parm,f_kvec,f_mat,title
REAL(DP)                                 :: parms(50),kx,ky,kz,norm,rgx
LOGICAL                                  :: lmat,ldos,lverb,lkvec,lmesh,ltest
!
INTEGER, PARAMETER                       :: ikv=11,ipar=9,ibnd=10, &
                                            ima=12,itst=15,ipri=8
REAL(DP)                                 :: dummyn     
!
parms  =0.0_DP
uplo   ='U'
jobz   ='N'
f_prot ='     '
f_band ='     '
f_parm ='     '
f_kvec ='     '
f_mat  ='     '
!
!***********************************************************************
!     Input of control data, prepare calculation
!***********************************************************************
!
CALL input_data
!
!***********************************************************************
!     Write protocol of control data
!***********************************************************************
!
OPEN(ipri,FILE=f_prot,FORM='formatted',STATUS='unknown' )
!
CALL protocol
!
!***********************************************************************
!     Calculate the eigenvalues
!***********************************************************************
!
!-----------------------------------------------------------------------
!     Allocation of arrays
!-----------------------------------------------------------------------
!
lwork=3*hdim
ALLOCATE(ham(hdim,hdim))
ALLOCATE(w(hdim))
ALLOCATE(work(3*hdim))
ALLOCATE(rwork(3*hdim-2))
!
IF(ldos) THEN
   ALLOCATE(bands(nkp,hdim))
   ALLOCATE(weight(hdim))
ENDIF 
!
!*******************************************************************
! header band structure file
!*******************************************************************
!
IF(lkvec) THEN
   WRITE(ibnd,10) nkp, hdim
ELSEIF(lmesh.AND.(dos.EQ.' TETRA')) THEN   
   WRITE(ibnd,11) nkp,hdim,ist,nteil,rgx   
ENDIF
!
IF(lmat) THEN
  jobz='V'
ENDIF 
!
IF(lverb.AND.(iverb.EQ.2)) THEN
    IF(.NOT.ldos) THEN 
       WRITE(il,*)
       WRITE(il,*) 'Calculated eigenvalues'
       WRITE(il,*) '----------------------'
    ELSE
       WRITE(il,*)
       WRITE(il,*)  ' List of bands'
       WRITE(il,*)  '-----------------' 
    ENDIF   
ENDIF  
!
!  loop k-vectors
!  
DO i=1,nkp   
   kx=kmesh(1,i)
   ky=kmesh(2,i)
   kz=kmesh(3,i)
!
   IF(lmat) WRITE(ima,10) i
! 
   CALL hamilton(ham,kx,ky,kz,parms)
   
! 
   CALL zheev(jobz,uplo,hdim,ham,hdim,w,work,lwork,rwork,info)
!
   WRITE(ibnd,75)  i, (w(k),k=1,hdim)
   IF(lverb.AND.(iverb.EQ.2)) WRITE(il,70) i, (w(k),k=1,hdim)
!
! calculation of weight factors for partial DOS
!      
   IF(ldos.AND.lmat) THEN   
      DO k=1,hdim      
         DO j=1,hdim
            weight(j)=ham(j,k)*CONJG(ham(j,k))
         ENDDO
         WRITE(ima,80) (weight(m),m=1,hdim)
      ENDDO      
   ENDIF              
ENDDO   
!
   CLOSE(ibnd)
   CLOSE(ima)
   CLOSE(il)
!
10 FORMAT(5I5)
11 FORMAT(4I5,F10.4)
20 FORMAT(8(F10.5,1x))
72 FORMAT(3x,I5,' | ',3(F10.4,2x))
70 FORMAT(3x,I5,' | ',8(F10.4,2x)/(11x,(8F10.4,2x)))
75 FORMAT(3x,I5,30(E14.6,2x))
80 FORMAT(30(E14.6,2x))
!
!***********************************************************************
!     Internal subroutines
!***********************************************************************
!       
      CONTAINS
!-----------------------------------------------------------------------
!     Input
!-----------------------------------------------------------------------
! 
       SUBROUTINE input_data
       IMPLICIT NONE
       INTEGER             :: cnt,io,npr,dosm,ios      
!     
       REAL(DP)            :: dummyn               
       CHARACTER(LEN=10)   :: t_code,dummy
!     
       OPEN(il,FILE='verbose.bnd',FORM='formatted',STATUS='unknown' )
       t_code='INIT  CARD'
       cnt=0
       ldos =.FALSE.
       lmat =.FALSE.
       lverb=.FALSE.
       lkvec=.FALSE.
       lmesh=.FALSE.
       ltest=.TRUE.
!
       DO WHILE (t_code.NE.'FINAL_CARD')
          READ(5,1,ADVANCE='NO') t_code
          WRITE(*,*) t_code
          If(lverb) THEN
             WRITE(il,*) t_code      
          ENDIF  
	      SELECT CASE(t_code)
	        CASE('.         ')
	           READ(5,5) dummy
!
!  start
!	           
	        CASE('PROBLEMTXT')
               READ(5,2) problem
            CASE('   VERBOSE')
               READ(5,*,IOSTAT=ios) iverb
               IF(ios.NE.0) THEN
                 WRITE(*,*) ' define verbosity level !'
                 STOP
               ELSE  
                 WRITE(*,*) ' verbosity level : ',iverb
                 lverb=.TRUE.
               ENDIF     
!	        
!  File names
!	        
            CASE('  F_MATELM')
	           READ(5,2) f_mat
               lmat=.TRUE.
               f_mat=ADJUSTL(f_mat)
               OPEN(ima	,FILE=f_mat,FORM='formatted',STATUS='unknown' )	     
            CASE('F_PROTOCOL')
	           READ(5,2) f_prot
               f_prot=ADJUSTL(f_prot)	
            CASE('F_BNDSTRUC')
	           READ(5,2) f_band
               f_band=ADJUSTL(f_band)
               OPEN(ibnd,FILE=f_band,FORM='formatted',STATUS='unknown' )	
            CASE('F_ETB_PARM')
	           READ(5,2) f_parm
               cnt=cnt+1
               f_parm=ADJUSTL(f_parm)
               OPEN(ipar,FILE=f_parm,FORM='formatted',STATUS='unknown' )
               nparm=-2
               DO
                 READ(ipar,*,IOSTAT=io) dummy,dummyn 
                 nparm=nparm+1
                 IF(io.LT.0) EXIT
               ENDDO 
               REWIND(ipar)                          
               READ(ipar,2) title
               DO i=1,nparm
                  READ(ipar,*) nam(i),parms(i)
               END DO   
               CLOSE(ipar)
            CASE('F_KVECTORS')             
               IF(lmesh) THEN
                  WRITE(ip,*) '   Conflict in k-mesh generation'
                  STOP
               ENDIF    
               lkvec=.TRUE.
	           READ(5,2) f_kvec
               f_kvec=ADJUSTL(f_kvec)
               OPEN(ikv,FILE=f_kvec,FORM='formatted',STATUS='unknown')
               READ(ikv,*) nkp
               ALLOCATE(kmesh(3,nkp))   
               DO i=1,nkp
                  READ(ikv,*)    k, kmesh(1,i), kmesh(2,i),kmesh(3,i)                 
               ENDDO 
               CLOSE(ikv)
!
!  Other data
!   
            CASE('  DOS_CALC')
               READ(5,*) dos
               dos=ADJUSTR(dos)
               IF(.NOT.((dos.EQ.' GAUSS').OR.(dos.EQ.' TETRA'))) THEN
                 WRITE(*,*) ' Argument of DOS_CALC has to be GAUSS or TETRA !' 
                 STOP
               ENDIF  
	           ldos =.TRUE.
	           ltest=.FALSE.      
            CASE('HAM_DIMENS') 
               READ(5,*)  hdim  
!
! K-mesh for DOS calculation
!                                      
            CASE(' KMESH_GEN')  
               IF(ltest) THEN
                  WRITE(ip,*) '   Card DOS_CALC before KMESH_GEN if DOS calculation'
                  STOP
               ENDIF    
               If(lkvec) THEN
                  WRITE(ip,*) '   Conflict in k-mesh generation'
                  STOP
               ENDIF    
               IF(.NOT.ldos) THEN
                 WRITE(ip,*) '   Fix method of DOS calculation for corresponding k-vector generation.'
                 STOP
               ENDIF
               IF(ldos.AND.(dos.EQ.'  GAUSS')) THEN
                  WRITE(ip,*) '   No special method for GAUSS implemented read k-vectors from file'
                  WRITE(ip,*) '   or use TETRA'
                  STOP
               ENDIF                              
               lmesh=.TRUE.
!
! K-mesh for tetrahedron method sc, fcc, bcc
!              
              IF(dos.EQ.' TETRA') THEN
                  READ(5,*) ist,nteil,rgx
!
! calculate k-point numbers
!                 
                  IF(ist.EQ.1) THEN
                    nkp=nteil*(1+nteil)*(2+nteil)/6
                  ELSEIF(ist.EQ.2) THEN
                    nkp=9*nteil**2+1+2*nteil*(8*nteil**2+7)/3   
                  ELSEIF(ist.EQ.3) THEN
                    nkp=(nteil+1)*(nteil+2)*(2*nteil+3)/6
                  ENDIF    
                  ALLOCATE(kmesh(3,nkp))
                  ALLOCATE(we(nkp))
                  IF(ist.EQ.1) CALL q3sc(nteil,nkp,rgx,kmesh,we)  
                  IF(ist.EQ.2) CALL q3fcc(nteil,nkp,rgx,kmesh,we)     
                  IF(ist.EQ.3) CALL q3bcc(nteil,nkp,rgx,kmesh,we)   
              ENDIF                      
! 	           
 	        CASE('FINAL_CARD')
                 IF(iverb.GT.1) WRITE(ip ,*) ' Input completed.'
      	    CASE DEFAULT
	           WRITE(ip,*) ' Keyword ' ,t_code,'  undefined -> STOP ' 
	           STOP
          END SELECT
       END DO   
!
1      FORMAT(A10,1x)
2      FORMAT(A60)
3      FORMAT(3I5)
4      FORMAT(3F10.6)
5      FORMAT(A10) 
51     FORMAT(A5) 
7      FORMAT(A3)
!
       END SUBROUTINE input_data
!
!-----------------------------------------------------------------------
!     Protocol
!-----------------------------------------------------------------------
!
      SUBROUTINE protocol
      IMPLICIT NONE
      INTEGER   :: k1,k2
      CHARACTER(LEN=5) :: TLAT(15)
!  Definition from Skriver LMTO. Not all is used here.      
      DATA TLAT/'   sc','  fcc','  bcc','  hex','   st','  bct',' trig', &
                '   so',' baco','  bco','  fco','   sm','  bcm',' tric', &
                'sowbz'/
!
      WRITE(ipri,1)
1     FORMAT(20x,' ********************************************'/   &
             20x,' *       Calculation of energy bands        *'/   &
             20x,' *------------------------------------------*'/   & 
             20x,' *       semi-empirical model (ETBM)        *'/   &
             20x,' ********************************************'/)
      WRITE(ipri,2) problem,f_prot,f_band,f_mat,f_parm,f_kvec
2     FORMAT(' Problem              : ',A60//, &
             ' Protocol             : ',A60/, &
             ' Band structure       : ',A60/, &
             ' Matrix elements      : ',A60/, &
             ' ETBM parameters      : ',A60/, &
             ' k-vectors            : ',A60/  //) 
      WRITE(ipri,*)
      WRITE(ipri,3)  hdim,nkp,nparm,ldos,lmat,title
3      FORMAT(1x,'Dimension of Hamiltonian   hdim  = ',I5,/&
              1x,'Number of k-vektors        nkp   = ',I5,/&
              1x,'Number of ETB-parameters   nparm = ' I5,/&
              1x,'Energy values for DOS calcuation   ',L5,/&    
              1x,'Weights for DOS calcuation         ',L5,//&             
              1x,'Source of ETB-parameters           ',A60//)
      
      IF(lmesh.AND.(dos.EQ.'TETRA')) THEN
         WRITE(ipri,*) 'Tetrahedron method (Dresden)'
         WRITE(ipri,*) '****************************' 
         WRITE(ipri,41) ist,tlat(ist),nteil,nkp,rgx
41       FORMAT(1x,'Lattice:                     LAT = ',I3,' - ',A5,  /&
                1x,'Subdivision k-mesh            NT = ',I3,/&
                1x 'Number of k-points           NKP = ',I3,/&
                1x,'Scaling                      RGX = 'F10.5 //)                 
      ENDIF                   
      WRITE(ipri,*)
      WRITE(ipri,*) 'Set of TB parameters'
      WRITE(ipri,*) '--------------------'
      DO i=1,nparm
         WRITE(ipri,120) i,ADJUSTR(nam(i)), parms(i)
      ENDDO  
     IF(lverb.AND.(iverb.GE.1)) THEN
         WRITE(il,*)
         WRITE(il,*) 'Table of k-points '
         WRITE(il,*)  '-----------------'
         DO i=1,nkp
            WRITE(il,70) i,kmesh(i,1),kmesh(i,2),kmesh(i,3)
         ENDDO   
      ENDIF
      WRITE(ipri,*)       
!
120   FORMAT(3x,I3,' | ',A25,F10.4)
70    FORMAT(3x,I5,' | ',12(F10.4,2x))
75    FORMAT(3x,'Number of PDOS : ',I3)
76    FORMAT(3x,'From orbital ',I3,' to orbital ',I3)
      CLOSE(ipri)
     END SUBROUTINE protocol 
!
END PROGRAM tb_bands
!
!***
