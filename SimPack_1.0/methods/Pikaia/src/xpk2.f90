
!****g* src/xpk2
! 
!  NAME
!   xpk2
!  AUTHOR
!   W. Hergert (modification)
!  TYPE
!   Test program - generalized least squares
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
      program xpk2

      implicit  none
      common/data/ xdat(200),ydat(200),scal(10),ndata
      real        :: xdat,ydat,scal,pi
      integer     :: ndata, n, seed, i, status,nfou,nparm
      real        :: ctrl(12), ctrl0(12),scal0(10),x(5), fit, fit2
      logical     :: l_final,l_data,l_scale
      external    :: fit2
      CHARACTER(LEN=60):: problem
      data pi/3.1415926536/
!
!     Set control variables (50 individuals for 200 generations
!     under Select-random-delete-worst reproduction plan)
!
      do 10 i=1,12
         ctrl0(i) = -1
   10 continue
      ctrl0(1)=50
      ctrl0(2)=200
      ctrl0(9)=0.0
      ctrl0(10)=3
      ctrl=ctrl0
! 
!     call input routine to overwrite ctrl 
!
      CALL ga_input(ctrl,problem,seed,nparm,nfou,l_final,l_data,l_scale)
!
!     read data
!      
      IF(l_data) THEN
         read(12,*) ndata
         DO i=1,ndata 
            read(12,*) xdat(i),ydat(i)
         ENDDO
      ELSE
         WRITE(*,*) ' no data provided -> STOP'
         STOP
      ENDIF        
!       
!    set the original  scale factors or use input file
!      
      scal0   =2.0
      scal0(5)=pi
      IF(l_scale) THEN
         DO i=1,nparm
            read(14	,*) scal(i)
         ENDDO
      ELSE
         scal=scal0
      ENDIF  
      WRITE(6,*) scal
!     Now call pikaia
!
      call rninit(seed)
      call pikaia(fit2,nparm,ctrl,x,fit,status) 
      !
!     rescale parameters
!
      WRITE(6,*) x
      x(1)=scal(1)*x(1)
      x(2)=scal(2)*x(2)
      x(3)=scal(3)*x(3)
      x(4)=scal(4)*x(4)
      x(5)=scal(5)*x(5)
!      
!
!     Print the results
!
      write(*,5) problem
 5    FORMAT(3x,A60)
      write(*,*) 
      write(*,*) '   status: ',status
      write(*,*) 
      write(*,*) '    final fit'
      write(*,25) 'x0,y0 : ', x(1),x(2)
      write(*,25) 'a,b   : ', x(3),x(4)
      write(*,25) 'theta : ', x(5)
25    FORMAT(3x,A10,3(E14.6,2x))      
      write(*,25) '   fit: ', fit 
      write(*,*)  '      ctrl:'
      write(*,20) ctrl
!
!     Export protocol to a file
!  
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*) '! least squares non-linear (see 5.1 p 35 ff. ) !'
       WRITE(8,*) '!----------------------------------------------!'
       WRITE(8,*)
       WRITE(8,5) problem
       WRITE(8,*)
       WRITE(8,*) '   vector of control data  (standard test) :'
       WRITE(8,*)
       WRITE(8,20) ctrl0
       WRITE(8,*)
       WRITE(8,*) '   vector of control data  ( actually used):'
       WRITE(8,*)
       WRITE(8,20) ctrl   
 20    FORMAT((3x,6(F12.5,2x)))
       WRITE(8,*)
       WRITE(8,*) '   status of calculation : ', status
       WRITE(8,*)
       WRITE(8,25) '   fit: ', fit 
       WRITE(8,*)
       WRITE(8,25) 'x0,y0 : ', x(1),x(2)
       WRITE(8,25) 'a,b   : ', x(3),x(4)
       WRITE(8,25) 'theta : ', x(5)         
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
       END
!
      function fatan(yy,xx)
!===========================================================
!     Returns arctangent in full circle
!===========================================================
      real yy,xx,fatan,a1,pi
      data pi/3.1415926536/

      a1=atan(yy/xx)
      if(xx.lt.0.) a1=a1+pi
      if(a1.lt.0.) a1=2.*pi+a1
      fatan=a1

      return
      end


