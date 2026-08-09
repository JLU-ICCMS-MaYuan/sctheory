!****m* src/dos_3d.f90
!
! NAME
!  dos_3d
! 
! AUTHOR
!  W. Hergert
!
! HISTORY
!  * 28.03.2022 checked
!  * 10.05.2022 a version that works finished
!  * 08.09.2022 verbosity levels checked
!  * 19.09.2022 implementation of tetrahedron method (Dresden)
!  * 15.10.2020 the implementation of tetrahedron method from 
!               Skrivers LMTO code not succesful, removed
!
! USAGE 
!  Calculation of the density of states and partial densities of 
!  of states from a FORTRAN band structure calculation by means of
!  an Hamiltonian exported from GTPack.
!
!
!
! DESCRIPTION
!  This program is part of SIMPack, a collection of simple FORTRAN 
!  routines to be used in connection with the Mathematica group theory 
!  package GTPack.
!
!  In this programme two versions of DOS calculations are realized:
!  (1)  calculation with Gaussian smearing
!  (2)  tetrahedron method
!        The tetrahedron method is extracted from Skrivers LMTO code)
!  
! INPUT
!  A set of tokens is defined to structure the input (CHARACTER(LEN=10)
!
!  File names :
!
!   'F_PROTOCOL'  -  filename protocol of the calculation
!   'F_BNDSTRUC'  -  filename bandstructure data    (in)
!   '  F_MATELM'  -  filname weight factors         (in)
!   ' F_DENSITY'  -  filename density of states     (out)
!   'F_PROTOCOL'  -  filename protocol calculation  (out)
!
!  Other control data  :
!
!   '   VERBOSE'  -  additional information , verbosity level
!   'PROBLMTEXT'    - headline for the calculation
!
!   'DOS_METHOD'  -  method DOS calculation (GAUSS, TETRA)
!   '  SMEARING'  -  half-width of Gaussian if GAUSS
!   'ENERG_AXIS'  -  definition of energy axis
!   'N_ELECTRON'  -  electron number
!   '  NORM_FAC'  -  normalisation factor
!   'BANDS_NUMB'  -  numbers of bands in calculation 
!   '  PART_DOS'  -  definition of the partial DOS
!   'FERI_ENERG'  -  Fermi energy (if known)
!   'FINAL_CARD'  -  last card of the input
!
! OUTPUT
!   The output can be found in files. Only a few information go to
!   the screen to control the correctness of the calculation.
!
! LITERATURE
!
!   TETRA
!   Lehmann, G. and Taut, M., On the Numerical Calculation of the 
!   Density of States and Related Properties,
!   physica status solidi (b)}, 54, 469-477 (1972)
!
! SOURCE
!
!
PROGRAM dos_3d
USE typedef
USE kpoints_3D
!
      IMPLICIT NONE
!
      REAL(DP), ALLOCATABLE, DIMENSION (:,:)     :: bands,dos,kmesh
      REAL(DP), ALLOCATABLE, DIMENSION (:,:,:)   :: weight,tl
      REAL(DP), ALLOCATABLE, DIMENSION (:)       :: en,dd,di,temp,we,a, &
                                                    e
      
      INTEGER(I4B)        :: iba,iw,ik,ie,j,k,nc,iverb
      INTEGER(I4B)        :: nkp,hdim,nenin,nb,nbu,nbo
      INTEGER(I4B)        :: kk,npart,ix,ipd,i,nteil,lat
      INTEGER(I4B)        :: npt(20),np(20,20)
! 
! npt and np control the partial DOS. The number of partial DOS is 
! limited to 20  
!    
      CHARACTER(LEN=60)   :: f_band,f_prot,f_dos,f_mat,problem
      CHARACTER(LEN=5)    :: dosm,npn(10)
      CHARACTER(LEN=6)    :: meth
      LOGICAL             :: lmat,ldos,lverb,lgauss,lpart
      
      REAL(DP)            :: emin,de,emax,ef,ga,zvbz,fac,fnorm,ns, &
                             nel(50),tmp,rgx,tlt
!   
      INTEGER, PARAMETER  :: ipr=8,idos=9,ibnd=10,imat=11, &
                             ndim=5000
      REAL(DP), PARAMETER :: splt=1.0D-4                                                               
!
!***********************************************************************
!     Input of control data,prepare calculation
!***********************************************************************
!
      CALL input_data  
!
!-----------------------------------------------------------------------
!     Read energy bands and weights
!-----------------------------------------------------------------------
!
      If(meth.EQ.' GAUSS') THEN
         READ(ibnd,*) nkp,hdim
      ELSEIF(meth.EQ.' TETRA') THEN
         READ(ibnd,*) nkp,hdim,lat,nteil,rgx
      ENDIF
!      
      ALLOCATE(bands(nkp,hdim))
!      
IF(lverb.AND.(iverb.GE.1)) THEN
       WRITE(il,*)
       WRITE(il,*)  ' List of bands'
       WRITE(il,*)  '-----------------'     
ENDIF  
      DO ik=1,nkp
         READ(ibnd,*) kk,(bands(ik,iba),iba=1,hdim)
         IF(lverb.AND.(iverb.EQ.2)) WRITE(il,70) kk,(bands(ik,iba),iba=1,hdim)
      ENDDO
!     
!
!       
      ALLOCATE(weight(nkp,hdim,hdim))  
      IF(lmat) THEN
        IF(lverb.AND.(iverb.EQ.2)) THEN
           WRITE(il,*)
           WRITE(il,*)  ' Weight factors'
           WRITE(il,*)  '-----------------'     
        ENDIF  
        DO ik=1,nkp
           READ(imat,*) kk
           IF(lverb.AND.(iverb.EQ.2)) WRITE(il,21) kk
           DO iba=1,hdim
              READ(imat,*) (weight(ik,iba,iw),iw=1,hdim)
              IF(lverb.AND.(iverb.EQ.2)) WRITE(il,71) (weight(ik,iba,iw),iw=1,hdim)
           ENDDO
        ENDDO
      ENDIF       
!
!-----------------------------------------------------------------------
!     Set energy axis
!-----------------------------------------------------------------------
!
      nenin=(emax-emin)/de
      ALLOCATE(en(nenin))
      DO ie=1,nenin
         en(ie)=emin+(ie-1)*de
      ENDDO
!
!***********************************************************************
!     Calculate DOS
!***********************************************************************
!
 WRITE(*,*) 'Allocate ',npart,nenin
      ALLOCATE(dd(nenin))
      ALLOCATE(di(nenin))      
      ALLOCATE(dos(nenin,npart+3))
      ALLOCATE(temp(nenin))
      dos=0.D0  
      di =0.D0 
      dd =0.D0 
!
IF(meth.EQ.' GAUSS') THEN
!
!-----------------------------------------------------------------------
!     Gaussian smearing
!-----------------------------------------------------------------------    
!
!  loop k-points and bands
!     
      DO ik=1,nkp
         DO iba=1,hdim
!           
!  total DOS
!
            temp=0.D0                          
            DO ie=1,nenin
                temp(ie)=temp(ie)+fnorm*        &
                EXP(-(en(ie)-bands(ik,iba))**2  &
                /2./ga/ga)/SQRT(2.0*pi)/ga
            ENDDO 
            DO ie=1,nenin
               dos(ie,1)=dos(ie,1)+temp(ie)
            ENDDO                  
!
!  all partial DOSes     
! 
            IF(lpart) THEN        
               DO ipd=1,npart
!
!   calculate t_l character  
!            
                  tlt=0.D0  
                  DO i=1,npt(ipd)
                     ix=np(ipd,i)
                     tlt=tlt+weight(ik,iba,ix)            
                  ENDDO        
                  temp=0.D0                          
                  DO ie=1,nenin
                     temp(ie)=temp(ie)+fnorm*tlt*      &
                     EXP(-(en(ie) -bands(ik,iba))**2   &
                     /2./ga/ga)/SQRT(2.0*pi)/ga
                  ENDDO  
!                  
!  sum up to partial DOSes and total DOS for control   
!               
                  DO ie=1,nenin
                     dos(ie,ipd+3)=dos(ie,ipd+3)+temp(ie)
                     dos(ie,3)=dos(ie,3)+temp(ie)
                  ENDDO
               ENDDO                 
            ENDIF         
!
!   end do loops k-vectors and bands
!                  
         ENDDO
      ENDDO   
!
!  normalize all results
!      
      DO ipd=1,npart+3
         DO ie=1,nenin
            dos(ie,ipd)=dos(ie,ipd)/nkp 
         ENDDO
      ENDDO         
!
ELSEIF(meth.EQ.' TETRA') THEN  
!

!
!-----------------------------------------------------------------------
!     Tetrahedron method (Dresden)
!-----------------------------------------------------------------------
! 
! allocate arrays
!
      ALLOCATE(kmesh(3,nkp))
      ALLOCATE(we(nkp))
      ALLOCATE(a(nkp))
      ALLOCATE(e(nkp))     
!
! generate kmesh
!       
      IF(lat.EQ.1) CALL q3sc (nteil,nkp,rgx,kmesh,we)  
      IF(lat.EQ.2) CALL q3fcc(nteil,nkp,rgx,kmesh,we)     
      IF(lat.EQ.3) CALL q3bcc(nteil,nkp,rgx,kmesh,we)  
!
! loop over bands
!       
      DO iba=nbu,nbo
         IF(lverb.AND.(iverb.EQ.2)) THEN
            WRITE(*,*) ' Band number    NB =  ',iba
         ENDIF    
!
!  calculate the total DOS
!  -----------------------
!
        temp=0.0
        DO k=1,nkp
           a(k)=1.0
           e(k)=bands(k,iba)
        ENDDO
!        
! DOS per band
!
        IF(lat.EQ.1) CALL GLBI(emin,de,emax,ef,splt,fnorm,ndim,nkp, &
                        nteil,nteil,e,a,e,kmesh,temp,tetkub)   
        IF(lat.EQ.2) CALL GLBI(emin,de,emax,ef,splt,fnorm,ndim,nkp, &
                        nteil,nteil,e,a,e,kmesh,temp,gltf)    
        IF(lat.EQ.3) CALL GLBI(emin,de,emax,ef,splt,fnorm,ndim,nkp, &
                        nteil,nteil,e,a,e,kmesh,temp,gltb)           
!
!      Summation  
!
       DO ie=1,nenin
          dos(ie,1)=dos(ie,1)+temp(ie)
       ENDDO
    ENDDO
!
!  calculate partial DOSes
!  -----------------------
!     
    IF(lpart) THEN                    
          DO iba=1,hdim              
!
!  all partial DOSes     
!     
               DO ipd=1,npart
!
!  calculate t_l character  
!            
                  a=0.0
                  DO ik=1,nkp
                     tlt=0.D0  
                     DO i=1,npt(ipd)
                        ix=np(ipd,i)
                        tlt=tlt+weight(ik,iba,ix)            
                     ENDDO
                     a(ik)=tlt
                     e(ik)=bands(ik,iba)
                  ENDDO           
                  temp=0.D0                          
!        
!  DOS per band
!
        IF(lat.EQ.1) CALL GLBI(emin,de,emax,ef,splt,fnorm,ndim,nkp, &
                        nteil,nteil,e,a,e,kmesh,temp,tetkub)   
        IF(lat.EQ.2) CALL GLBI(emin,de,emax,ef,splt,fnorm,ndim,nkp, &
                        nteil,nteil,e,a,e,kmesh,temp,gltf)    
        IF(lat.EQ.3) CALL GLBI(emin,de,emax,ef,splt,fnorm,ndim,nkp, &
                        nteil,nteil,e,a,e,kmesh,temp,gltb)  
!                                 
!  sum up to partial DOSes and total DOS for control    
!              
                  DO ie=1,nenin
                     dos(ie,ipd+3)=dos(ie,ipd+3)+temp(ie)
                     dos(ie,3)=dos(ie,3)+temp(ie)
                  ENDDO
              ENDDO
            ENDDO
                              
      ENDIF         
!   
ENDIF
!
!-----------------------------------------------------------------------
!  calculate integrated DOS and Fermi energy, same for both methods
!-----------------------------------------------------------------------                     
!
      DO ie=1,nenin
         dd(ie)=dos(ie,1)
      ENDDO
      CALL dqsf(de,dd,di,nenin)
      DO ie=1,nenin
         dos(ie,2)=di(ie)
      ENDDO 
!
      k=0
      DO ie=1,nenin
         IF((di(ie)-zvbz).LT.0.D0) THEN
            k=ie 
          ELSE
            EXIT
         ENDIF   
         Write(*,*) k,en(ie), (di(ie)-zvbz) 
      END DO
      WRITE(*,*) k
      ef=emin+(k-1)*de+(zvbz-di(k))/(di(k+1)-di(k))*de
 123 CONTINUE     
!
!-----------------------------------------------------------------------
!     Output in file
!-----------------------------------------------------------------------
!
      WRITE(idos,22) emin,de,emax,ef,zvbz
      IF(lpart) THEN
         nc=npart+3
         WRITE(idos,21) nc,nenin
         DO iba=1,nc
            WRITE(idos,21) iba
            WRITE(idos,22) (dos(ie,iba),ie=1,nenin)
         END DO  
      ELSE
         nc=2
         WRITE(idos,21) nc,nenin
         DO iba=1,2
            WRITE(idos,21) iba
            WRITE(idos,22) (dos(ie,iba),ie=1,nenin)
         END DO  
      ENDIF      
!
      CLOSE(idos)      
!
!***********************************************************************
!     Protocol of input data
!***********************************************************************
!
      CALL protocol
!
!***********************************************************************
!
8     FORMAT(1X,'FERMI-ENERGIE  Ef=',E14.6/)
21    FORMAT(16I5)
22    FORMAT(8(F10.5,1X))
26    FORMAT(15E13.6)
70    FORMAT(3x,I5,' | ',8(F10.4,2x)/(11x,(8F10.4,2x)))
71    FORMAT(3x,8(F10.4,2x)/(11x,(8F10.4,2x)))
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
       INTEGER                     :: lattv,pnt,ip,ix,iy,iz, &
                                      i,j,ndist,ios
       CHARACTER*10                :: t_code,dummy
       CHARACTER*3                 :: lattice
       DOUBLE PRECISION, PARAMETER :: EPS=1.D-05
!
!      verbose information to file
!
       OPEN(il,FILE='verbose.dos',FORM='formatted',STATUS='unknown' )
       !
       t_code='INIT  CARD'
       f_mat ='     '
       lmat  =.FALSE.
       ldos  =.TRUE.
       lverb =.FALSE.
       lpart =.FALSE.
       npart = 0
!
       DO WHILE (t_code.NE.'FINAL_CARD')
          READ(5,1,ADVANCE='NO') t_code
          If(lverb) THEN
             WRITE(il,*) t_code
          ENDIF   
! tokens to screen          
	      SELECT CASE(t_code)
	        CASE('.         ')
	           READ(5,5) dummy
!
! Start
!
            CASE('   VERBOSE')
               READ(5,*,IOSTAT=ios) iverb
               IF(ios.NE.0) THEN
                 WRITE(*,*) ' define verbosity level !'
                 STOP
               ELSE  
                 WRITE(*,*) ' verbosity level : ',iverb
                 lverb=.TRUE.
               ENDIF        
            CASE('PROBLEMTXT')
               READ(5,2) problem
!
! File names
!	           
            CASE('F_PROTOCOL')
	           READ(5,2) f_prot
               f_prot=ADJUSTL(f_prot)
               OPEN(ipr,FILE=f_prot,FORM='formatted',STATUS='unknown' )
            CASE('F_BNDSTRUC')
	           READ(5,2) f_band
               f_band=ADJUSTL(f_band)
               OPEN(ibnd,FILE=f_band,FORM='formatted',STATUS='unknown' )
            CASE(' F_DENSITY')
	           READ(5,2) f_dos
	           ldos=.TRUE.
               f_dos=ADJUSTL(f_dos)
               OPEN(idos,FILE=f_dos,FORM='formatted',STATUS='unknown' )
            CASE('  F_MATELM')
	           READ(5,2) f_mat
               lmat=.TRUE.
               f_mat=ADJUSTL(f_mat)
               OPEN(imat,FILE=f_mat,FORM='formatted',STATUS='unknown' )
!
! Other control data
! 
            CASE('ENERG_AXIS')
	           READ(5,*) emin,de,emax
            CASE('N_ELECTRON')
	           READ(5,4) zvbz         
            CASE('BANDS_NUMB')
	           READ(5,*) nb,nbu,nbo
            CASE(' FERMI_ENE')
	           READ(5,*) ef
            CASE('  NORM_FAC')
	           READ(5,*) fnorm   
            CASE('DOS_METHOD')
               READ(5,*) meth
	           meth=ADJUSTR(meth)  
	           IF(meth.EQ.' GAUSS') ldos=.FALSE.	            
	           IF(meth.EQ.' TETRA') ldos=.FALSE.  
	           IF(ldos) THEN
	              WRITE(il,*) 'Wrong method DOS calculation! -> STOP'
	              WRITE(*,*)  'Wrong method DOS calculation! -> STOP'
	              STOP
	           ENDIF
	        CASE('  PART_DOS')
               READ(5,*) npart
	           lpart=.TRUE.
! 
! Read information about partial DOS: name, which orbitals to sum	 
!          
	           DO i=1,npart
	              READ(5,*) npn(i),npt(i)
	              READ(5,*) (np(i,j),j=1,npt(i))
	           ENDDO
            CASE('  SMEARING')
               lgauss=.TRUE.
	           READ(5,4) ga
 	        CASE('FINAL_CARD')
               IF(iverb.GT.1) WRITE(ipr,*) ' Input competed! '
      	     CASE DEFAULT
	            WRITE(6,*) ' Keyword ' ,t_code,'  undefined -> STOP ' 
	            STOP
          END SELECT
       END DO   
!
!***********************************************************************
!
1      FORMAT(A10,1x)
2      FORMAT(A60)
3      FORMAT(3I5)
4      FORMAT(3F10.6)
5      FORMAT(A10) 
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
      INTEGER  :: ii,j
      WRITE(ipr,1)
1     FORMAT(20x,' ********************************************'/   &
             20x,' *  Calculation of density of states        *'/   &
             20x,' *------------------------------------------*'/   & 
             20x,' *       semi-empirical model (ETBM)        *'/   &
             20x,' ********************************************'/)
      WRITE(ipr,"(1x,'Problem : '//A80/)") problem
      WRITE(ipr,2) f_prot,f_band,f_mat,f_dos
2     FORMAT(' Protocol             : ',A60/, &
             ' Band structure       : ',A60/, &
             ' Weight factors       : ',A60/, &
             ' Density of States    : ',A60//)
      WRITE(ipr,3)  emin,de,emax,zvbz,ef,fnorm,nkp,hdim,nb,nbu,nbo
3     FORMAT(1x,'Energy axis             emin,de,emax = ',(3F10.6,1x),/&
             1x,'Number of electrons             zvbz = ',F10.6,/ &
             1x,'Fermi energy                      ef = ',F10.6,/ &
             1x,'Normalization factor           fnorm = ',F10.6,/ &
             1x,'Number of k-points               nkp = ',I5,//   &
             1x,'Dimension Hamiltonian           hdim = ',I5,/    &
             1x,'Total number of bands             nb = ',I5,/    &
             1x,'Number of lowest band u.c.       nbu = ',I5,/    &
             1x,'Number of last band u.c.         nbo = ',I5,//)
!
      IF(lgauss) THEN
            WRITE(ipr,4) ga
4           FORMAT(1x,'width for Gauss smearing          ga = ',F10.6//)
      ELSE
            WRITE(ipr,5) nteil,lat
5           FORMAT(1x,'Subdivision k-mesh                   = ',I3,/&
                   1x,'Lattice type                     lat = ',I3,//)  
      ENDIF 
!        
      IF(lpart) THEN
           WRITE(ipr,*) 'Calculation of partial DOS :'
           WRITE(6,*) npart
           DO ii=1,npart
              WRITE(6,*) npn(ii),npt(ii)
              WRITE(ipr,6)  npn(ii),(np(ii,j),j=1,npt(ii))
           ENDDO       
      ENDIF
6     FORMAT(1x,A5,2x,10(I3,2x))   
!   
      END SUBROUTINE protocol 
!
END PROGRAM dos_3d

