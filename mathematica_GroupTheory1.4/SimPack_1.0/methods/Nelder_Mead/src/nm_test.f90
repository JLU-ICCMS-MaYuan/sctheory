!****m* src/nm_test.f90
!
! NAME
!  nm_test.f90
!
! AUTHOR
!  W. Hergert 
!
! USE
! Test of Nelder-Mead routine asa047.f90
!
! DESCRIPTION
!
!  The_original test program asa047_test.f90 was given by John Burkardt.
!  In this version the tests can be modified individually. The input data
!  are coming from files. Standard values correspondto the original test.
!  The structure of the program is taken from the TB fitting procedure.
!
! INPUT
!  A set of tokens is defined to structure the input (CHARACTER(LEN=10)
!
!   'PROBLMTEXT'    - headline for the calculation
!   '   F_START'  -  filname approximation start    (in)
!   '   F_FINAL'  -  filename approximation final   (out)
!   'F_PROTOCOL'  -  filename protocol calculation  (out)
!
!   '   VERBOSE'  -  level of verbosity
!   'FINAL_CARD'  -  last card of the input
!
!***********************************************************************
! 
! SIMPLEX METHOD (MINIM):
! ----------------------
!
!    '  SIM_STEP' - initial stepsizes
!    ' SIM_MAXFN' - maximum nnumber of function evaluations
!    'SIM_IPRINT' - print control
!                    <0 no printing
!                    =0 parameters and function values after convergence
!                    >0 like =0, but with progress report
!    '  SIM_STOP' - stopping criterion
!    ' SIM_NLOOP' - Stopping rule applied after NLOOP function 
!                   evaluations (standard: 2*NOP)
!    ' SIM_IQUAD' - =1 if fitting for quadratic surfce (recommended)
!                   =0 if not
!    '  SIM_SIMP' - criterion for expanding the simplex
! 
!
!***********************************************************************
!    
! SOURCE
!
PROGRAM nm_test
!
USE typedef
!
CHARACTER(LEN=60)   :: f_prot,f_start,f_final,title,problem
CHARACTER(LEN=8)    :: meth
LOGICAL             :: lros,lpow,lquart,lheli,lin,lverb,lfint
INTEGER(I4B)        :: iverb
INTEGER, PARAMETER  :: ipri=8,ist=9,ifin=10,ibnd=12,icnt=15
!
CHARACTER(LEN= 8),ALLOCATABLE,DIMENSION(:)  :: names
INTEGER(I4B)                                :: nparm
INTEGER(I4B)                                :: i
!
! Simplex Nelder-Mead (MINIM)
!
REAL(DP), EXTERNAL                          :: roesenbrock, quaric, helical, powell
REAL(DP), ALLOCATABLE,DIMENSION(:)          :: start0,start,step0,step,xmin
REAL(DP)                                    :: reqm0,reqm,diff1,diff2
INTEGER(I4B)                                :: conv,conv0,cont,cont0,icount,numres,ifault
!  
CALL timestamp()                              
CALL input_data
!
ALLOCATE(start0(nparm))
ALLOCATE(start(nparm))
ALLOCATE(step0(nparm))
ALLOCATE(step(nparm))
ALLOCATE(xmin(nparm))
ALLOCATE(names(nparm))
!
! READ start names, parameters and steps
!
If(lin) THEN
   DO i=1,nparm
      READ(ist,*) names(i), start(i), step(i) 
   ENDDO
   start0=start
ENDIF    
!
! Perform the tests
!  
IF(lros) THEN
!
! Rosenbrock
!
! original test data
!
  start0(1:nparm)= (/-1.2D+00,1.0D00/)
  step0(1:nparm) = (/1.0D+00,1.0D+00/)
  names(1:nparm) =(/'P(1)','P(2)' /)
  reqm0         = 1.0D-08
  conv0          = 10
  cont0          = 500
  IF(.NOT.lin) THEN
    start=start0
    step =step0
    reqm =reqm0 
    conv  =conv0
    cont  =cont0
  ENDIF
!    
  diff1 = rosenbrock ( start0 )
  call nelmin ( rosenbrock, nparm, start, xmin, diff2, reqm, step, &
    conv,cont, icount, numres, ifault )
!
ELSEIF(lpow) THEN
!
! Powell
!
  start0(1:nparm) = (/  3.0D+00, - 1.0D+00,   0.0D+00,   1.0D+00 /)
  step0(1:nparm)  = (/ 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /)
  names(1:nparm) =(/'P(1)','P(2)','P(3)','P(4)' /)
  reqm0 = 1.0D-08
  conv0 = 10
  cont0 = 500
  IF(.NOT.lin) THEN
    start=start0
    step =step0
    reqm =reqm0 
    conv  =conv0
    cont  =cont0
  ENDIF
!    
  diff1 = powell ( start0 )
  call nelmin ( powell, nparm, start, xmin, diff2, reqm, step, &
    conv, cont, icount, numres, ifault )
!
ELSEIF(lquart) THEN
!
! Quartic
!
  start0(1:nparm) = 1.0D+00
  step0(1:nparm)  = 1.0D+00
  names(1:nparm) =(/' P(1)',' P(2)',' P(3)',' P(4)',' P(5)', &
                    ' P(6)',' P(7)',' P(8)',' P(9)','P(10)'/)
  reqm0 = 1.0D-08
  conv0 = 10
  cont0 = 2000
  IF(.NOT.lin) THEN
    start=start0
    step =step0
    reqm =reqm0 
    conv  =conv0
    cont  =cont0
  ENDIF
!
  diff1= quartic ( start0 )
  call nelmin ( quartic, nparm, start, xmin, diff2, reqm, step, &
    conv, cont, icount, numres, ifault )
!
ELSEIF(lheli) THEN
!
! Helical
!
  start0(1:nparm) = (/ - 1.0D+00,   0.0D+00,   0.0D+00 /)
  step0(1:nparm)  = (/ 1.0D+00, 1.0D+00, 1.0D+00 /)
  names(1:nparm) =(/'P(1)','P(2)','P(3)'/)
  reqm0 = 1.0D-08
  conv0 = 10
  cont0 = 500
  IF(.NOT.lin) THEN
    start=start0
    step =step0
    reqm =reqm0 
    conv  =conv0
    cont  =cont0
  ENDIF
!    
  diff1 = helical ( start0 )
  call nelmin ( helical, nparm, start, xmin, diff2, reqm, step, &
    conv,cont, icount, numres, ifault )
!
ENDIF
! 
!
!  final result of the fit
!
!   
IF(lfint) THEN
   DO i=1,nparm
     WRITE(ifin,*) xmin(i)
   ENDDO
ENDIF
!
CALL timestamp
CALL protocol

!
!**********************************************************************
!
!     internal subroutines
!
!**********************************************************************
!
      CONTAINS
!
!  input of data
!
       SUBROUTINE input_data
       IMPLICIT NONE
       CHARACTER*10                :: t_code,dummy
       CHARACTER(LEN=8)            :: meth
       INTEGER(I4B)                :: cnt
!
!      verbose information to file
!
       OPEN(il,FILE='verbose.fit',FORM='formatted',STATUS='unknown' )
       !
       t_code='INIT  CARD'
       cnt=0
       lros    =.FALSE.
       lpow    =.FALSE.
       lquart  =.FALSE.
       lheli   =.FALSE.
       lin     =.FALSE.
       lfint   =.FALSE.
       f_start ='---'
       f_final ='---'
!
       DO WHILE (t_code.NE.'FINAL_CARD')
          READ(5,1,ADVANCE='NO') t_code
          If(lverb) THEN
             WRITE(il,*) t_code
          ELSE
             WRITE(*,*) t_code
          ENDIF   
! tokens to screen          
	      SELECT CASE(t_code)
	        CASE('.         ')
	           READ(5,3) dummy
!
!  file names and problem  
!	           
	        CASE('PROBLEMTXT')
               READ(5,2) problem
               cnt=cnt+1
            CASE('F_PROTOCOL')
	           READ(5,2) f_prot
               cnt=cnt+1
               f_prot=ADJUSTL(f_prot)          
            CASE('   F_START')
                lin=.TRUE.
	           READ(5,2) f_start
               cnt=cnt+1
               f_start=ADJUSTL(f_start)
               OPEN(ist,FILE=f_start,FORM='formatted',STATUS='unknown' )
            CASE('   F_FINAL')
	           READ(5,2) f_final
	           lfint=.TRUE.
               cnt=cnt+1
               f_final=ADJUSTL(f_final)
               OPEN(ifin,FILE=f_final,FORM='formatted',STATUS='unknown' )
!
!  general
!            

            CASE('   VERBOSE')
               READ(5,*) iverb
               lverb=.TRUE.
               cnt=cnt+1       
            CASE('     INPUT')
               READ(5,*) dummy
               lin=.TRUE.
               cnt=cnt+1                  
            CASE(' PARM_NUMB')
	           READ(5,*) nparm
               cnt=cnt+1
            CASE('  FIT_FUNC')
               READ(5,4) meth
	           meth=ADJUSTR(meth)
	           IF(meth.EQ.'ROSENBRO') THEN
	             lros=.TRUE.
	           ELSEIF(meth.EQ.'  POWELL') THEN 
	             lpow=.TRUE.
	           ELSEIF(meth.EQ.' QUARTIC') THEN
	             lquart=.TRUE.
	           ELSEIF(meth.EQ.' HELICAL') THEN
	             lheli=.TRUE.
	           ENDIF    
!             
! control for simplex method
!
            
            CASE(' SIM_MAXFN')
                cnt=cnt+1
                READ(5,*) cont
            CASE('  SIM_STOP')
                READ(5,*) reqm
                cnt=cnt+1
            CASE(' SIM_NLOOP')
                cnt=cnt+1
                READ(5,*) conv          
!
 	        CASE('FINAL_CARD')
               IF(cnt.GT.120) THEN
                  WRITE(il,*) ' Information in input missing  -> STOP '
                  WRITE(  6,*) ' Information in input missing  -> STOP '
                  STOP
               ENDIF  
      	     CASE DEFAULT
      	        WRITE(il,*) ' Keyword ' ,t_code,'  undefined -> STOP ' 
	            WRITE(6,*) ' Keyword ' ,t_code,'  undefined -> STOP ' 
	            STOP
          END SELECT
       END DO   
!
1      FORMAT(A10,1x)
2      FORMAT(A60)
3      FORMAT(A10)
4      FORMAT(A8)
!
       END SUBROUTINE input_data
!
!
! protocol
!
      SUBROUTINE protocol
      IMPLICIT NONE
      INTEGER :: i,l
      OPEN(ipri,FILE=f_prot,FORM='formatted',STATUS='unknown' )
      WRITE(ipri,1)
1     FORMAT(3x,'******************************************'/ &
             3x,'*  Tests of simplex method according to  *'/ &
             3x,'*              Nelder - Mead             *'/ &
             3x,'******************************************') 
      WRITE(ipri,"(1x,'Problem : '//A80/)") problem
!
! file names
!      
      WRITE(ipri,2) f_prot,f_start,f_final
2     FORMAT('   Protocol              : ',A60/, &
             '   parameters start      : ',A60/, &
             '   parameters final      : ',A60//)
!
! general
!
      WRITE(ipri,*)
      WRITE(ipri,3) nparm
3     FORMAT(3x,'Number of parameters  : ',I3)
      WRITE(ipri,*)              
!
! simplex Nelder-Mead (nelim)
!
      
      WRITE(ipri,*) ' Test is :'
      IF(lros)   WRITE(ipri,*) '          -> Rosenbrock'
      IF(lpow)   WRITE(ipri,*) '          -> Powell'
      IF(lquart) WRITE(ipri,*) '          -> Quartic'
      IF(lheli)  WRITE(ipri,*) '          -> Helical'
      WRITE(ipri,*)
      WRITE(ipri,4) diff1,diff2
4     FORMAT(3x,'Initial deviation      : ',E12.4/, & 
             3x,'Final deviation        : ',E12.4//) 
      WRITE(ipri,*) 

      WRITE(ipri,92) cont,reqm,conv,icount,numres,ifault
92    FORMAT(3x,'               In/OUT '/, &
             3x,'--------------------------------------'/, &
             3x,' MAXFN  ',I5/, &     
             3x,'  STOP  ',E12.4/, &   
             3x,' NLOOP  ',I5/, & 
             3x,'icount  ',I5 /, &
             3x,'numres  ',I5 /, &  
             3x,'ifault  ',I5)
      WRITE(ipri,*)
!
      WRITE(ipri,*)      
      WRITE(ipri,7)
7     FORMAT(3x,'********************************************'/ &
             3x,'*         Result of the fit                *'/ &
             3x,'*     parameter    initial      final      *'/ &        
             3x,'********************************************') 
      WRITE(ipri,*)
      DO l=1,nparm
           WRITE(ipri,8) l,ADJUSTR(names(l)),start0(l),xmin(l)
      ENDDO
8     FORMAT(3x,I2,2x,A12,2x,F10.5,2x,E14.6)   

      END SUBROUTINE protocol 

END PROGRAM nm_test

