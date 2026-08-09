SUBROUTINE ga_input(ctrl,problem,seed,nparm,nfou,l_final,l_data,l_scale)
!**********************************************************************
!
!     internal subroutines
!
!**********************************************************************
!
       IMPLICIT NONE
       REAL                        :: ctrl(12),xmin,xmax
       CHARACTER*10                :: t_code,dummy
       CHARACTER(LEN=60)           :: f_prot,f_start,f_final,problem,f_data,f_scale
       INTEGER                     :: cnt,ir,seed,nparm,nfou
       LOGICAL                     :: l_final,l_data,l_scale
!
       l_final=.FALSE.
       l_data =.FALSE.
       l_scale=.FALSE.	
!
!      verbose information to file
!
!
       t_code='INIT  CARD'
       cnt=0
       ir=5
       write(*,*) 'drin'
!    
       DO WHILE (t_code.NE.'FINAL_CARD')
          READ(5,1,ADVANCE='NO') t_code
          WRITE(*,*) t_code
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
               OPEN(8,FILE=f_prot,FORM='formatted',STATUS='unknown' ) 
            CASE('    F_DATA')
	           READ(5,2) f_data
	           l_data=.TRUE.
               cnt=cnt+1
               f_data=ADJUSTL(f_data)
               OPEN(12,FILE=f_data,FORM='formatted',STATUS='unknown' )   
            CASE('   F_SCALE')
	           READ(5,2) f_scale
	           l_scale=.TRUE.
               cnt=cnt+1
               f_scale=ADJUSTL(f_scale)
               OPEN(14,FILE=f_scale,FORM='formatted',STATUS='unknown' )    
            CASE('   F_FINAL')
               l_final=.TRUE.
	           READ(5,2) f_final
               cnt=cnt+1
               f_final=ADJUSTL(f_final)
               OPEN(10,FILE=f_final,FORM='formatted',STATUS='unknown' )
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
            CASE('  GEN_SEED')	       
               cnt=cnt+1
               READ(ir,*) seed
            CASE(' PARM_NUMB')	       
               cnt=cnt+1
               READ(ir,*) nparm
            CASE('   FOURIER')	       
               cnt=cnt+1
               READ(ir,*) nfou           
!
 	        CASE('FINAL_CARD')
               IF(cnt.GT.120) THEN
                  WRITE(  6,*) ' Information in input missing  -> STOP '
                  STOP
               ENDIF  
      	     CASE DEFAULT
	            WRITE(6,*) ' Keyword ' ,t_code,'  undefined -> STOP ' 
	            STOP
          END SELECT
       END DO   
1      FORMAT(A10,1x)
2      FORMAT(A60)
3      FORMAT(A10)
4      FORMAT(A8)       
END SUBROUTINE ga_input
