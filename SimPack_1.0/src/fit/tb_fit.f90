!****f* src/tb_fit.f90
!
! NAME
!  tb_fit.f90
!
! AUTHOR
!  W. Hergert
!
! USE
!  This is the driver routine for estimation of TB paramter sets
!  from ab initio band structures.
! 
!  The following methods are implemented:
!   * least squares fit           (DNLQ1)
!   * genetic algorithm           (Picaia)
!   * simplex method Nelder-Mead  (asa047)
! 
!   It is possible to comine the methods, i.e. to takethe output of one method as 
!   input of another method.
!  arbitrary initial values. It can be used after a genetic algorithm
!  or a simplex search.
!
! INPUT
!  A set of tokens is defined to structure the input (CHARACTER(LEN=10)
!
!   'PROBLMTEXT'  -  headline for the calculation
!   'F_ABINITIO'  -  filename data for fitting      (in)
!   '   F_START'  -  filname approximation start    (in)
!   '   F_FINAL'  -  filename approximation final   (out)
!   'F_PROTOCOL'  -  filename protocol calculation  (out)
!   ' F_FIXPARM'  - information on fixed parameters (in) 
!
!   '   VERBOSE'  -  additional information  
!   'FIT_METHOD'  -  method used in the fit
!   ' DATA_FORM'  -  structue of ab initio data
!
!   'HAM_DIMENS'  -  dimension of the Hamiltonian
!   '   VERBOSE'  -  level of verbosity
!   'FINAL_CARD'  -  last card of the input
!   'WEIGHT_BND'  -  filename for weights of bands in fitting
!
!***********************************************************************
!
! LEAST SQARES (DNLQ1) :
! ----------------------
!
!   ' DNLQ1_EPS'  - stopping criterion relative change of parameters
!   ' DNLQ1_EFX'  - stopping criterion defec
!   ' DNLQ1_EGR'  - stopping criterion norm of gradient
!   ' DNLQ1_ITM'  - number of interations
!   ' DNLQ1_IPR   - print control
!                   = 0 no printing
!                   = 1 print start and final values
!                   = 2 as 1, but additional each iteration step
!                   = 3 as 2, but additional all intermediate values
!                   = 4 as 3, but additional new Papproximation vector
!   ' DNLQ1_IST'  - control value
!                   = 0 store duplicate Jacobi matrix
!                   = 1 recalculate Jacobi matrix in each step
!
!***********************************************************************
!
! GENETIC ALGORITHM (PICAIA) :
! ----------------------------
! 
!   ' GEN_INDIV'  - number of individuals in population (standard 100)
!   ' GEN_GENER'  - number of generations (standard 500)
!   '  GEN_NUMB'  - number of significant digits, i.e.  number ogf genes
!   ' GEN_CROSS'  - crossover probability rate (default 0.85)
!   'GEN_MUTMOD'  - mutation mode 1/2=steady/variable (default is 2)
!   'GEN_MUTINI'  - initial mutation rate (default 0.005)
!   'GEN_MUTMIN'  - minimum mutation rate => 0.0 (default 0.0005) 
!   'GEN_MUTMAX'  - maximum mutation rate <= 1.0 (default 0.25) 
!   'GEN_RELFIT'  - relativ fintness diferential (default is 1.)
!   'GEN_REPROD'  - reproduction plan 
!                   1 - fllgenerational replacement
!                   2 - steady-state-replace-random
!                   3 - steady state replace-worst (default)
!   'GEN_ELITIS'  - elitism flag  0/1 = off/on (default 0)  
!   ' GEN_PRINT'  - printed output 0/1/2=none/minimal/verbose (default 0)
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
! DESCRIPTION
!    
! SOURCE
!
MODULE common_fit
REAL, SAVE    :: f(200),t(200),sigma
INTEGER,SAVE :: ndata
END MODULE common_fit
!
PROGRAM tb_fit
!
USE typedef
USE tbdata
!
! COMMON Block to fix parameters of the set
!
REAL(DP)           :: vfix
INTEGER(I4B)       :: ifp,nump,nfix
COMMON /pfix/ vfix(40),ifp(40),nump,nfix,lfix
!
! general variables
!
CHARACTER(LEN=60)   :: f_abinit,f_prot,f_start,f_final,title,problem, &
                       f_fitc, f_fixp,dummyt,fixtitle,f_simplex,      &
                       f_scale,simptitle
CHARACTER(LEN=6)    :: meth
CHARACTER(LEN= 8)   :: code
LOGICAL             :: lsimp,lverb,lgauss,lgen,lparm,lwgt,lfint,lcont, &
                       lfix,lstep
INTEGER, PARAMETER  :: ipri=8,ist=9,ifin=10,ibnd=12,icnt=15,icnp=16,&
                       isim=25,isc=28
!
CHARACTER(LEN= 8),ALLOCATABLE,DIMENSION(:)  :: names
INTEGER(I4B)                                :: hamd,ios,iverb,nact
INTEGER(I4B)                                :: i,isel

! Least squares (DNLQ1)
!
REAL(DP)                                    :: diff1,diff2,epsrel,efx,egrad, &
                                               epsrel0,efx0,egrad0
INTEGER(I4B)                                :: itm,iprint,ie,itm0,iprint0,ie0
EXTERNAL                                    :: F1,tb_ls
REAL(DP), ALLOCATABLE,DIMENSION(:)          :: parms,x,fx,dwork,parms0,parms1
!
! Simplex Nelder-Mead (MINIM)
!
INTEGER, PARAMETER                          :: rk=kind(1.0D+00)
REAL(DP), EXTERNAL                          :: tb_simplex
REAL(DP), ALLOCATABLE,DIMENSION(:)          :: step,xn
REAL(DP)                                    :: simp,stopcr,stepi
INTEGER                                     :: maxfn,nloop,iquad,ifault, numres,icount
LOGICAL                                     :: first
!
! genetic algorithm (PIKAIA)
!
REAL                                        :: ctrl(12),ctrli(12),fit
REAL, ALLOCATABLE,DIMENSION(:)              :: xg
INTEGER                                     :: status
! 
!**********************************************************************
!
! Define standard values for the methods and read the starting 
! approximation
!
!**********************************************************************
!
!
! standard values gentic algorithm
!                                  
DATA ctrl/100.,500.,6.,0.85,2.0,0.005,0.0005,0.25,1.,3.,0.,0./
!
! 
! standard values simplex method
!   
   maxfn =10000
   nloop =100
   stopcr=1.0D-10
! 
!  Input of control data
!  ---------------------
!
      CALL input_data      
! 
!
!  Read initial parameter set
!  --------------------------
!
      ios   =0
      nparm =-2
      DO WHILE (ios.EQ.0)
         nparm=nparm+1	
         READ(ist,*,IOSTAT=ios) dummyt  
      ENDDO
      REWIND(ist)
!
! Allocate arrays 
!     parms0 - parameters from file
!     parms1 - final result
!     names  - names of parameters
!     step   - necessary for simplex method
!
      ALLOCATE(parms(nparm))
      ALLOCATE(parms0(nparm))
      ALLOCATE(parms1(nparm))
      ALLOCATE(names(nparm))
!
      READ(ist,11) title
11    FORMAT(A60)
      DO i=1,nparm
         READ(ist,*) names(i),parms(i) 
      ENDDO  
      CLOSE(ist)
      parms0=parms
! 
!**********************************************************************
!
! Prepare for restricted parameter sets and weights
!
!**********************************************************************
!
! Check if some parameters are fixed
! ----------------------------------
!
      IF(lfix) THEN
        ios   =0
        nfix  =-1
        DO WHILE (ios.EQ.0)
           nfix=nfix+1
           READ(icnp,*,IOSTAT=ios) dummyt  
        ENDDO
        REWIND(icnp)
        IF(nfix.GT.40) THEN
          WRITE(*,*) 'ERROR : number of fixed parameters exceeds 40 -> STOP'
          STOP
        ENDIF 
        READ(icnp,*) fixtitle
        nfix=nfix-1
        DO i=1,nfix
           READ(icnp,*) ifp(i),vfix(i)
        ENDDO  
        CLOSE(icnp)
!
!  Check if fixed to input parameters
!  (number in ifp is negative)
!       
        DO i=1,nfix
           IF(ifp(i).LT.O) THEN
             ifp(i)=ABS(ifp(i)) 
             vfix(i)=parms(ifp(i))
           ENDIF  
        ENDDO
!   
!  copy fixed parameters in parms0, i.e.
!  parms0 reflects the starting approximation wiith fixed parameters
!        
        l=1
        DO i=1,nparm
           IF(i.EQ.ifp(l)) THEN
              parms0(i)=vfix(l)
              l=l+1
           ENDIF
        ENDDO    
        nump=nparm
        nact=nparm-nfix
        ALLOCATE(x(nact))
        isel=-1
        CALL parmfix(x,parms,isel)
!
! nothing fixed
!        
      ELSE
        nact=nparm
        nump=nparm
        nfix=0  
        ALLOCATE(x(nact))
        x=parms
      ENDIF    
!      
!   weights for bands
!          
        IF(.NOT.lwgt) THEN
          ALLOCATE(weight(hdim))
          weight=1.0
        ENDIF            
!

      IF(iverb.GE.2) THEN
         WRITE(*,*)
         WRITE(*,*) '   Number of parameters         : ',nump
         WRITE(*,*) '   Number of fixed parameters   : ',nfix
         WRITE(*,*) '   Acutual number of parameters : ',nact
         WRITE(*,*)
      ENDIF   
!
      IF(iverb.GE.2) THEN
         WRITE(*,*) '   Initial approximation'
         DO i=1,nact
            WRITE(*,12) i,x(i)
         ENDDO
      ENDIF
12    FORMAT(3x,I3,2x,F10.5)     
      WRITE(*,*)
      WRITE(*,*) 
!
!**********************************************************************
!
!  Read energies for the fit
!
!**********************************************************************
! 
      CALL abinitio(f_abinit,code,iverb)
!
!**********************************************************************
!
! Fit with the different methods
!
!**********************************************************************
! 
!   least squares (DNLQ1)
!   ---------------------
!
      IF(lgauss) THEN
!
! form the vector for the deviations
!
        neq=SUM(sel)
        ALLOCATE(fx(neq))
!
! deviation with initial  parameter set
!  
        CALL tb_ls(nact,x,neq,fx,F1,ier)
!
        diff1=DOT_PRODUCT(fx,fx)/SQRT(neq*1.0_DP)
        IF(lverb.AND.iverb.EQ.2) THEN
          WRITE(il,*)
          WRITE(il,*) 'initial difference : ',diff1
          WRITE(*,*)  'initial difference : ',diff1         
        ENDIF  
!
        ALLOCATE(dwork(3*nact+2*neq+2*nact*neq))
!
        CALL DNLQ1(nact,x,neq,tb_ls,f1,epsrel,efx,egrad,itm,iprint,ie,dwork)
!
! deviation with final parameter set
! 
        CALL tb_ls(nact,x,neq,fx,F1,ier)
        diff2=DOT_PRODUCT(fx,fx)/SQRT(neq*1.0_DP) 
!
        IF(lverb.AND.iverb.EQ.2) THEN
          WRITE(il,*) 'final difference   : ',diff2
          WRITE(*,*)  'final difference   : ',diff2
        ENDIF          
!              
      ELSEIF(lgen) THEN
!
!  genetic algorithm
!  -----------------
!       
        ctrli=ctrl
        ALLOCATE(xg(nparm))
        xg=0.

! Rescale the parameters to interval 0 <= x_i <= 1
!
        IF(lscal) THEN
          ALLOCATE(scal(nparm))
          DO i=1,nparm
             READ(isc,*) k,scal(i)
          ENDDO   
          DO i=1,nparm
             xg(i)=scal(i)*parms(i)
          ENDDO 
        ELSE        
          DO i=1,nparm
             xg(i)=(parms(i)-xmin)/(xmax-xmin)
          ENDDO   
        ENDIF
!            
        IF(lverb.AND.(iverb.EQ.2)) THEN
          WRITE(il,*)
          WRITE(il,*) 'Scaling of parameters'
          DO i=1,nparm
             WRITE(il,62) i,ADJUSTR(names(i)),parms(i),xg(i)     
          ENDDO
          WRITE(il,*)
        ENDIF          
62    FORMAT(3x,I2,2x,A10,2x,F12.5,2x,F12.5)            
!          
        diff1= 1./tb_ga(nparm,xg)
!
        call rninit(seed)
        CALL pikaia(tb_ga,nparm,ctrl,xg,fit,status)
!
! deviation with final  parameter set
!
	    diff2=1./tb_ga(nparm,xg)
!
! Rescale finally parameter set
!
        DO i=1,nparm
           IF(lscal) THEN
              xg(i)=xg(i)/scal(i)
           ELSE   
              xg(i)=xg(i)*(xmax-xmin)+xmin
           ENDIF             
        ENDDO 

        x=xg
!
      ELSEIF(lsimp) THEN
!
!  simplex method
!  --------------
! 
        ALLOCATE(step(nparm))
        IF(lstep) THEN
          READ(isim,*) simptitle
          DO i=1,nparm
             READ(isim,*)  k,step(i)
          ENDDO
        ELSE  
          step=1.0
        ENDIF   
!  
! deviation with initial  parameter set
!
       ALLOCATE(xn(nparm))
!
       diff1= tb_simplex(parms)
! 
       CALL nelmin(tb_simplex,nparm,parms,xn,diff2,stopcr,step,nloop, &
                      maxfn,icount,numres,ifault)
!
       x=xn
     ENDIF       
!   
!**********************************************************************
!
! Output of results
! -----------------
!
! full parameter set after fitting
!
  IF(lfix) THEN
    CALL parmfix(x,parms1,+1)
  ELSE
    parms1=x  
  ENDIF  
!
OPEN(ipri,FILE=f_prot,FORM='formatted',STATUS='unknown' )
CALL protocol
!
!  final result of the fit
!   
IF(lfint) THEN
   DO i=1,nparm
     WRITE(ifin,*) parms1(i)
   ENDDO
   CLOSE(ifin)
ENDIF
!
IF(lcont) THEN
   WRITE(icnt,*) title
   DO i=1,nparm
     WRITE(icnt,101) names(i),parms1(i)
   ENDDO
   CLOSE(icnt)
ENDIF
101 FORMAT(A10,2x,E14.6)  
!

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
       lgauss =.FALSE.
       lgen   =.FALSE.
       lsimp  =.FALSE.
       lverb  =.FALSE.
       lparm  =.FALSE.
       lwgt   =.FALSE.
       lfint  =.FALSE.
       lcont  =.FALSE.
       lfix   =.FALSE.
       lstep  =.FALSE.
       lscal  =.FALSE.
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
            CASE('F_ABINITIO')
	           READ(5,2) f_abinit
               cnt=cnt+1
               f_abinit=ADJUSTL(f_abinit)              
            CASE('   F_START')
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
            CASE(' F_FITCONT')
	           READ(5,2) f_fitc
               cnt=cnt+1
               lcont=.TRUE.
               f_fitc=ADJUSTL(f_fitc)
               OPEN(icnt,FILE=f_fitc,FORM='formatted',STATUS='unknown' )   
            CASE(' F_FIXPARM')
	           READ(5,2) f_fixp
               cnt=cnt+1
               lfix=.TRUE.
               f_fixp=ADJUSTL(f_fixp)
               OPEN(icnp,FILE=f_fixp,FORM='formatted',STATUS='unknown' )   
!             
!  general
!            

            CASE('   VERBOSE')
               READ(5,*) iverb
               lverb=.TRUE.
               cnt=cnt+1                  
            CASE('HAM_DIMENS')
	           READ(5,*) hdim
               cnt=cnt+1
            CASE('FIT_METHOD')
               READ(5,4) meth
	           meth=ADJUSTR(meth)
	           IF(meth.EQ.' LEASTSQ') THEN
	             lgauss=.TRUE.
	           ELSEIF(meth.EQ.'  GENALG') THEN 
	             lgen=.TRUE.
	           ELSEIF(meth.EQ.' SIMPLEX') THEN
	             lsimp=.TRUE.
	           ENDIF    
	           cnt=cnt+1   
	        CASE(' DATA_FORM')
               READ(5,4) code
               code=ADJUSTL(code)
	           cnt=cnt+1   
	        CASE('WEIGHT_BND')
	           lwgt=.TRUE.
	           cnt=cnt+1
	           ALLOCATE(weight(hdim))
	           READ(5,*) (weight(i),i=1,hdim)
!
!  control for least squares
! 
            CASE(' DNLQ1_EPS')
               cnt=cnt+1
	           READ(5,*) epsrel
	           epsrel0=epsrel	        
	        CASE(' DNLQ1_EGR')	      
               cnt=cnt+1
	           READ(5,*) egrad
	           egrad0=egrad
	        CASE(' DNLQ1_EFX')
               cnt=cnt+1
	           READ(5,*) efx
	           efx0=efx
	        CASE(' DNLQ1_ITM')
               cnt=cnt+1
	           READ(5,*) itm
	           itm0=itm
	        CASE(' DNLQ1_IPR')
               cnt=cnt+1
	           READ(5,*) iprint
	           iprint0=iprint
	        CASE(' DNLQ1_IST')	       
               cnt=cnt+1
	           READ(5,*) ie
	           ie0=ie	       
!
! control for genetic algorithm
! default values are set, input used to override default
!
            CASE(' GEN_INDIV')	       
               cnt=cnt+1
               READ(ir,*) ctrl(1)
! CASE(' GEN_INDIV')	
            CASE(' GEN_GENER')       
               cnt=cnt+1
               READ(ir,*) ctrl(2)
            CASE('  GEN_NUMB')	       
               cnt=cnt+1
               READ(ir,*) ctrl(3)
            CASE(' GEN_CROSS')	       
               cnt=cnt+1
               READ(ir,*) ctrl(4)
            CASE('GEN_MUTMOD')	       
               cnt=cnt+1
               READ(ir,*) ctrl(5)
            CASE('GEN_MUTINI')       
               cnt=cnt+1
               READ(ir,*) ctrl(6)
            CASE('GEN_MUTMIN')	       
               cnt=cnt+1
               READ(ir,*) ctrl(7)
            CASE('GEN_MUTMAX')	       
               cnt=cnt+1
               READ(ir,*) ctrl(8)
            CASE('GEN_RELFIT')	       
               cnt=cnt+1
               READ(ir,*) ctrl(9)
            CASE('GEN_REPROD')	       
               cnt=cnt+1
               READ(ir,*) ctrl(10)
            CASE('GEN_ELITIS')	       
               cnt=cnt+1
               READ(ir,*) ctrl(11)
            CASE(' GEN_PRINT')	       
               cnt=cnt+1
               READ(ir,*) ctrl(12)
            CASE('  GEN_XMIN')	     
               cnt=cnt+1
               READ(ir,*) xmin
            CASE('  GEN_XMAX')	  
               cnt=cnt+1
               READ(ir,*) xmax
            CASE('  GEN_SEED')	       
               cnt=cnt+1
               READ(ir,*) seed
            CASE(' GEN_SCALE')	       
               cnt=cnt+1
               lscal =.TRUE.
               READ(5,2) f_scale
               f_scale=ADJUSTL(f_scale)
               OPEN(isc,FILE=f_scale,FORM='formatted',STATUS='unknown' ) 
             
! control for simplex method
!
            
            CASE(' SIM_MAXFN')
                cnt=cnt+1
                READ(5,*) maxfn
            CASE('  SIM_STOP')
                READ(ir,*) stopcr
                cnt=cnt+1
            CASE(' SIM_NLOOP')
                cnt=cnt+1
                READ(ir,*) nloop     
            CASE('  SIM_STEP')
                cnt=cnt+1
                lstep=.true.       
                READ(5,2) f_simplex
                f_simplex=ADJUSTL(f_simplex)
                OPEN(isim,FILE=f_simplex,FORM='formatted',STATUS='unknown' )       
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
! protocol
!
      SUBROUTINE protocol
      INTEGER :: l     
      WRITE(ipri,12)
12    FORMAT(3x,'******************************************'/ &
             3x,'*  Fit TB-model against ab initio data   *'/ &
             3x,'******************************************') 
      WRITE(ipri,*)
      WRITE(ipri,"(1x,'Problem : '//A80/)") problem
!
! file names
!      
      IF(.NOT.lcont) THEN
         f_fitc='     '
      ENDIF 
      IF(.NOT.lfix) THEN
         f_fixp='     '
      ENDIF   
      WRITE(ipri,22) f_prot,f_abinit,f_start,f_final,f_fitc,f_fixp
22    FORMAT('   Protocol                 : ',A60/, &
             '   ab initio data           : ',A60/, &
             '   TB parameters start      : ',A60/, &
             '   TB parameters final      : ',A60/, &
             '   TB parameters fit cont.  : ',A60/, &
             '   TB paprameters fixed     : ',A60//)
!
! general
!
      WRITE(ipri,*)
      WRITE(ipri,33) hdim,code
33    FORMAT(3x,'Dimension  Hamiltonian : ',I3/, &
             3x,'FORM of input data     : ',A8)
      WRITE(ipri,*)   
!
! parameters before fitting
!   
      WRITE(ipri,5)
5     FORMAT(3x,'Initial parameter set'/ &
             3x,'----------------------')                       
      WRITE(ipri,'(3x,A60)') ADJUSTL(title)
      WRITE(ipri,*)
      IF(lfix) THEN
         l=1        
         WRITE(ipri,*)
         DO i=1,nparm
            IF(i.EQ.ifp(l)) THEN
               WRITE(ipri,61) i,ADJUSTR(names(i)),parms(i),' parameter fixed to : ',vfix(l)            
               l=l+1
            ELSE 
               WRITE(ipri,6) i,ADJUSTR(names(i)),parms(i)
            ENDIF   
        ENDDO
      ELSE
        DO i=1,nparm
           WRITE(ipri,6) i,ADJUSTR(names(i)),parms(i)
        ENDDO
      ENDIF    
6     FORMAT(3x,I2,2x,A10,2x,F12.5,2x,F12.5)
61    FORMAT(3x,I2,2x,A10,2x,F12.5,A22,F12.5)    
      WRITE(ipri,*)
      WRITE(ipri,*)             
!
! least squares (DNLQ1)
!             
      IF(lgauss) THEN
         WRITE(ipri,*) '  Fit method             : least squares fit'
         WRITE(ipri,*)
         WRITE(ipri,*) '  Weight factors for bands :'  
         IF(lwgt) THEN           
            WRITE(ipri,41) (weight(i),i=1,hdim)
41          FORMAT(3x,(10(F8.2,2x)))
         ELSE
           WRITE(ipri,*) '  All factors are 1.0'  
         ENDIF    
         WRITE(ipri,*)
         WRITE(ipri,*)
         WRITE(ipri,4) diff1,diff2
4        FORMAT(3x,'Initial deviation      : ',E12.4/, & 
                3x,'Final deviation        : ',E12.4/) 
             
         WRITE(ipri,*)
         WRITE(ipri,*)
         WRITE(ipri,9) epsrel0,epsrel,efx0,efx,egrad0,egrad, &
                       itm0,itm,iprint0,iprint,ie0,ie
9        FORMAT(3x,'               In     |    Out '/, &
                3x,'--------------------------------------'/, &
                3x,'EPSREL  ',E12.4,'  |   ',E12.4/, &   
                3x,'   EFX  ',E12.4,'  |   ',E12.4/, &   
                3x,' EGRAD  ',E12.4,'  |   ',E12.4/, &   
                3x,' ITMAX  ',I5,'         |   ',I5/, &   
                3x,'IPRINT  ',I5,'         |   ',I5/, &   
                3x,'    IE  ',I5,'         |   ',I5)
                WRITE(ipri,*)
                WRITE(ipri,*)   
!
! simplex Nelder-Mead (nelim)
!
      ELSEIF(lsimp) THEN
        WRITE(ipri,*) '  Fit method             : simplex Nelder-Mead'
        WRITE(ipri,*)
        WRITE(ipri,4) diff1,diff2
        WRITE(ipri,*) 

        WRITE(ipri,92) maxfn,stopcr,nloop,icount,numres,ifault
92      FORMAT(3x,'               In/OUT '/, &
                3x,'--------------------------------------'/, &
                3x,' MAXFN  ',I5/, &     
                3x,'  STOP  ',E12.4/, &   
                3x,' NLOOP  ',I5/, & 
                3x,'icount  ',I5 /, &
                3x,'numres  ',I5 /, &  
                3x,'ifault  ',I5)
                WRITE(ipri,*)
                WRITE(ipri,*)
! 
! genetic algorithm (PIKAIA)  
!      
      ELSEIF(lgen) THEN
        WRITE(ipri,*) '  Fit method             : genetic algorithm'
        WRITE(ipri,*)
        WRITE(ipri,4) diff1,diff2
        WRITE(ipri,*)
        WRITE(ipri,*) ' final fitness : ',fit
        WRITE(ipri,*) '        status : ',status
        WRITE(ipri,*)
        WRITE(ipri,91) INT(ctrli(1)),INT(ctrl(1)), &
                       INT(ctrli(2)),INT(ctrl(2)), &
                       INT(ctrli(3)),INT(ctrl(3)), &
                       ctrli(4), ctrl(4), &
                       INT(ctrli(5)),INT(ctrl(5)), &
                       ctrli(6),ctrli(6),ctrli(7),ctrl(7), &
                       ctrli(8),ctrl(8),ctrli(9),ctrl(9), &
                       INT(ctrli(10)),INT(ctrl(10)), &
                       INT(ctrli(11)),INT(ctrl(11)), &
                       INT(ctrli(12)),INT(ctrl(12)),xmin,xmax,seed
91      FORMAT(3x,'               In     |    Out '/, &
               3x,'--------------------------------------'/, &
               3x,' GEN_INDIV',I5,'         |   ',I5/, &   
               3x,' GEN_GENER',I5,'         |   ',I5/, &   
               3x,'  GEN_NUMB',I5,'         |   ',I5/, &   
               3x,' GEN_CROSS',E12.4,'  |   ',E12.4/, &   
               3x,'GEN_MUTMOD',I5,'         |   ',I5/, & 
               3x,'GEN_MUTINI',E12.4,'  |   ',E12.4/, &   
               3x,'GEN_MUTMIN',E12.4,'  |   ',E12.4/, &  
               3x,'GEN_MUTMAX',E12.4,'  |   ',E12.4/, &   
               3x,'GEN_RELFIT',E12.4,'  |   ',E12.4/, &   
               3x,'GEN_REPROD',I5,'         |   ',I5/, &   
               3x,'GEN_ELITIS',I5,'         |   ',I5/, &   
               3x,' GEN_PRINT',I5,'         |   ',I5//, &
               3x ' GEN_XMIN ',E12.4/, &
               3x ' GEN_XMAX ',E12.4/, &
               3x ' GEN_SEED ',I8)
      WRITE(ipri,*)    
      WRITE(ipri,*)   
    ENDIF  
!   
      WRITE(ipri,7)
7     FORMAT(3x,'********************************************'/ &
             3x,'*         Result of the fit                *'/ &
             3x,'*     parameter    initial      final      *'/ &        
             3x,'********************************************') 
      WRITE(ipri,*)
      DO l=1,nparm
           WRITE(ipri,8) l,ADJUSTR(names(l)),parms0(l),parms1(l)
      ENDDO
8     FORMAT(3x,I2,2x,A12,2x,F10.5,2x,E14.6)   
      CLOSE(ipri)
      END SUBROUTINE protocol 
!
END PROGRAM tb_fit
!***

