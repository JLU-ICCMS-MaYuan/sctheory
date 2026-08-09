!****g* src/xpk1a
! 
!  NAME
!   xpk1a
!  AUTHOR
!   W. Hergert (modification)
!  TYPE
!   Test program - linear least square problem
!  PACKAGE
!   Pikaia
!  USAGE
!   Test of genetic algorithm
!  DESCRIPTION
!   * The originalcode from Pikaia distribution is used.
!   * An input subprogram is used with the tokens of the TB fit program.
!     This allows to modify easily the control parameters. 
!  MODIFICATION HISTORY
!   14th July 2022  : modification of the original test
!
! SOURCE
!  
      program xpk1a
      implicit none        
      common/data/ f(200),t(200),scal(100),sigma,ndata
      character(LEN=60)              :: problem
      real                           :: f,t,sigma,scal
      integer                        :: seed, i, status, nparm, ndata,nfou
      real                           :: ctrl(12), ctrl0(12),scal0(100),fit1a
      logical                        :: l_final,l_data,l_scale
      real, allocatable,dimension(:) :: x
      external                       :: fit1a  
!
!     Set control variables (evolve 50 individuals over 100
!     generations, use defaults values for other input parameters)
!
      scal0   =1.
      scal0(1)=10.0
      scal0(2)=100.0
      do i=1,12
         ctrl0(i) = -1
      enddo
      ctrl0(1)=50
      ctrl0(2)=100
      ctrl=ctrl0
! 
!     call input routine to overwrite ctrl and scl
!
      CALL ga_input(ctrl,problem,seed,nparm,nfou,l_final,l_data,l_scale)
!   
      IF(l_data) THEN
         read(12,*) ndata
         DO i=1,ndata 
            read(12,*) t(i),f(i)
         ENDDO
         sigma=5.
      ELSE
         WRITE(*,*) ' no data provided -> STOP'
         STOP
      ENDIF  
!
!     scaling
!       
      IF(l_scale) THEN
         DO i=1,nparm
            read(14	,*) scal(i)
         ENDDO
      ELSE
         scal=scal0
      ENDIF           
      ALLOCATE(x(nparm))
!
!     First, initialize the random-number generator
!
      call rninit(seed)
!
!     Now call pikaia
!
      call pikaia(fit1a,nparm,ctrl,x,f,status)
!
!     Print the results , Parameter are rescaled
!      
!     rescale parameters
!
      x(1)=scal(1)*x(1)
      x(2)=scal(2)*x(2)
!      
!
!     Print the results
!
      write(*,*) problem
      write(*,*)
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: '
      WRITE(*,10) f
10    FORMAT(5(F12.6,2x))      
      write(*,*) '   ctrl:'
      write(*,20) ctrl
!
!     Export protocol to a file
!  
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*) '!  least squares linear (see 5.1 p 35 ff. )    !'
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*)
       WRITE(8,21) problem
21     FORMAT(3x,A60)       
       WRITE(8,*) '   vector of control data  (standard test) :'
       WRITE(8,*)
       WRITE(8,20) ctrl0
       WRITE(8,*)
       WRITE(8,*) '   vector of control data  ( actually used):'
       WRITE(8,*)
       WRITE(8,20) ctrl
 20    FORMAT(3x,6f12.5/3x,6f12.5)
       WRITE(8,*)
       WRITE(8,*) '   status of calculation : ', status
       WRITE(8,*) '   result of the fit     : ', x
!
!      Output fit
!
       IF(l_final) THEN
          DO i=1,nparm
             WRITE(10,30) i, x(i)
          ENDDO
       ENDIF   
 30    FORMAT(I5,2x,F14.6)  
!              
      end 


