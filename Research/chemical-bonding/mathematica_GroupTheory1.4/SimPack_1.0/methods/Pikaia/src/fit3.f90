
      function fit3(n,x)
!=============================================================
!     Fitness function for circle fitting problem using 
!     a robust estimator based on the Hough transform
!     (Sect. 5.5)
!=============================================================
      implicit none
      common/data/ xdat(200),ydat(200),sigma,ndata
      integer      n,ndata,i
      real         x(n),fit3,xdat,ydat,sigma,x0,y0,distdat,  &
                   r1,r2,sum
!---------- 1. rescale input variables:
!           0 <= x(1,2) <= 3
      x0    = x(1)*3.
      y0    = x(2)*3.
!---------- 2. compute merit function
      r1=1.-sigma/2.
      r2=1.+sigma/2.
      sum=0.
      do i=1,ndata
!     (a) compute distance d_j from center to data point
         distdat=sqrt((x0-xdat(i))**2+(y0-ydat(i))**2)
!     (b) increment Hough merit function
         if(distdat.ge.r1.and.distdat.le.r2) sum=sum+1.
      enddo
!---------- 3. equate fitness to Hough merit function
      fit3=sum

      return
      end


