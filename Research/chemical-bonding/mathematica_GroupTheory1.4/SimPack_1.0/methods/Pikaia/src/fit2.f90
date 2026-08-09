
      function fit2(n,x)
!=============================================================
!     Fitness function for ellipse fitting problem (Sect. 5.4)
!
!     Ellipse parameters are:
!     x(1)=x_0, x(2)=y_0, x(3)=a, x(4)=b, x(5)=theta_0
!=============================================================
      implicit none
      common/data/ xdat(200),ydat(200),scal(10),ndata
      integer      n,ndata,i
      real         x(n),xdat,ydat,fit2,fatan,x0,y0,a2,b2, &
                   theta0,distdat,angdat,rthj,sum,pi,scal
      parameter    (pi=3.1415926536)
!---------- 1. rescale input variables:
!           0 <= x(1...4) <= 2, 0 <= x(5) <= pi
      x0    = x(1)*scal(1)
      y0    = x(2)*scal(2)
      a2    =(x(3)*scal(3))**2
      b2    =(x(4)*scal(4))**2
      theta0= x(5)*scal(5)
!---------- 2. compute merit function
      sum=0.
      do i=1,ndata
!     (a) compute distance d_j from center to data point
         distdat=sqrt((x0-xdat(i))**2+(y0-ydat(i))**2)
!     (b) compute angle theta_j of segment center---data point
         angdat=fatan((ydat(i)-y0),(xdat(i)-x0))-theta0
!     (c) compute radius r(\theta_j) of ellipse at that angle
         rthj=sqrt(a2*b2/(b2*cos(angdat)**2+a2*sin(angdat)**2))
!     (d) increment DR merit function
         sum=sum+(distdat-rthj)**2
      enddo
!---------- 3. equate fitness to inverse of DR merit function
      fit2=1./sum

      return
      end

