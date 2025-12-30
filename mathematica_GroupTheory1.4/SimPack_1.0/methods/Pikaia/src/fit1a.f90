
      function fit1a(n,x)
!========================================================
!     Fitness function for linear least squares problem
!     (Sect. 5.2)
!========================================================
      implicit none
      common/data/ f(200),t(200),scal(100),sigma,ndata
      integer      n,i,ndata
      real         x(n),fit1a,f,t,sigma,a,b,sum,scal
!---------- 1. rescale input variables:
      a=x(1)*scal(1)
      b=x(2)*scal(2)
!---------- 2. compute chi**2
      sum=0.
      do  i=1,ndata
         sum=sum+ ( (a*t(i)+b-f(i))/sigma )**2
      enddo
!---------- 3. define fitness
      fit1a=1./sum

      return
      end


