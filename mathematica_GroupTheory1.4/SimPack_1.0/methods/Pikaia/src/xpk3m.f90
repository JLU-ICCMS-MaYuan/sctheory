
      program xpk3
!============================================================
!     Driver program for circle fitting problem using
!     a robust estimator based on the Hough transform
!     (Sect. 5.5)
!============================================================
      implicit  none
      integer   n, seed, i, status
      parameter (n=2)
      real      ctrl(12), x(n), f, fit3
      external  fit3
!
!     First, initialize the random-number generator
!
      seed = 123456
      call rninit(seed)
!
!     Read in synthetic data
!
      call finit
!
!     Set control variables (use 
!     
      do i=1,12
         ctrl(i) = -1
      enddo

!     Now call pikaia
      call pikaia(fit3,n,ctrl,x,f,status)
!
!     rescaling
!
      x=3.*x
!     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      x: ',x
      write(*,*) '      f: ',f
      write(*,20) ctrl
 20   format(   '    ctrl: ',6f11.6/10x,6f11.6)

      end
!***************************************************************
      subroutine finit
!===============================================================
!     Reads in synthetic dataset
!===============================================================
      implicit     none
      common/data/ xdat(200),ydat(200),sigma,ndata
      real         xdat,ydat,sigma
      integer      ndata,i

      open(1,file='test3.dat',form='formatted')

      read(1,*) ndata
      DO i=1,ndata
         read(1,*) xdat(i),ydat(i)
      ENDDO
!
!     Same error estimate in x and y for all data points
      sigma=0.05

      return
      end
