!*********************************************************************
      function twod(n,x)
!=====================================================================
!     Compute sample fitness function (2-d landscape)
!=====================================================================
      implicit none
!
!     Input:
      integer n
      real x(n)
!
!     Output
      real twod
!
!     Constant
      real pi,sigma2
      integer nn
      parameter (pi=3.1415926536,sigma2=0.15,nn=9)
!
!     Local
      real rr

      if (x(1).gt.1..or.x(2).gt.1.) stop
      rr=sqrt( (0.5-x(1))**2+ (0.5-x(2))**2)
      twod=cos(rr*nn*pi)**2 *exp(-rr**2/sigma2)

      return
      end
