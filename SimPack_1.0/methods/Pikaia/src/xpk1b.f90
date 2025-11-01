!****g* src/xpk1b
! 
!  NAME
!   xpk1a
!  AUTHOR
!   W. Hergert (modification)
!  TYPE
!   Test program - non-linear least square problem
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
      program xpk1b
      implicit none        
      common/data/ f(200),t(200),scal(100),sigma,ndata,nfou
      character(LEN=60)              :: problem
      real                           :: f,t,sigma,scal
      integer                        :: seed, i, status, nparm, ndata,nfou
      real                           :: ctrl(12), ctrl0(12),scal0(100),fit1b
      real                           :: nyp,fit
      logical                        :: l_final,l_data,l_scale
      real, allocatable,dimension(:) :: x
      external                       :: fit1b 
      
    
!      real      ctrl(12), x(n), f, fit1b,nyp,sigma,t,f0
!
!     set the original control data 
      do i=1,12
         ctrl0(i) = -1
      enddo
      ctrl0(1)=50
      ctrl=ctrl0
! 
!     call input routine to overwrite ctrl 
!
      CALL ga_input(ctrl,problem,seed,nparm,nfou,l_final,l_data,l_scale)
!   
      IF(l_data) THEN
         read(12,*) ndata
         DO i=1,ndata 
            read(12,*) t(i),f(i)
         ENDDO
         sigma=5.
         nyp=2.0*(t(2)-t(1))
      ELSE
         WRITE(*,*) ' no data provided -> STOP'
         STOP
      ENDIF  
!       
!    set the original  scale factors or use input file
!      
      scal0   =1.0
      scal0(1)=10.0
      scal0(2)=100.0
      DO i=1,nfou
         scal0(3*i)=100.0
      ENDDO
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
      call pikaia(fit1b,nparm,ctrl,x,fit,status)
!
!     Print the results , Parameter are rescaled
!      
!     rescale parameters
!
      x(1)=scal(1)*x(1)
      x(2)=scal(2)*x(2)
      DO i=1,nfou
         x(3*i)  =x(3*i)*scal(3*i)
         x(3*i+1)=x(3*i+1)*(50.0-nyp)+nyp
      ENDDO   
!      
!
!     Print the results
!
      write(*,*) problem
      write(*,*) ' status: ',status
      write(*,*) '      x: '
      write(*,25) x(1),x(2)
      DO i=1,nfou
         write(*,25) x(3*i),x(3*i+1),x(3*i+2)
      ENDDO
25    FORMAT(3x,3(E14.6,2x))      
      write(*,*) '      f: '
      WRITE(*,10) fit
10    FORMAT(5(F12.6,2x))      
      write(*,*) '   ctrl:'
      write(*,20) ctrl
!
!     Export protocol to a file
!  
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*) '! least squares non-linear (see 5.1 p 35 ff. ) !'
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*)
       WRITE(8,*) problem
21     FORMAT(3x,A60)
       WRITE(8,*)       
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
       WRITE(8,*) '   final fitness         : ', fit
       WRITE(8,*) '   result of the fit     : ' 
       WRITE(8,*)
       WRITE(8,20) x(1),x(2)
       DO i=1,nfou
         WRITE(8,20) x(3*i),x(3*i+1),x(3*i+2)
       ENDDO
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
