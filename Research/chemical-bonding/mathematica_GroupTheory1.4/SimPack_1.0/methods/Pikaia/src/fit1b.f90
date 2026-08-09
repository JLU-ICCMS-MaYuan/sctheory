
      function fit1b(n,x)
!==========================================================
!     Fitness function for non-linear least squares problem
!     (Sect. 5.3)
!==========================================================
      implicit none
      common/data/ f(200),t(200),scal(100),sigma,ndata,nfou
      integer      n,nfou,mmax,ndata,i,j
      real         x(n),fit1b,f,t,sigma,amp,per,a,b, &
                   sum,sum2,nyp,pi,scal
      parameter    (mmax=10,pi=3.1516926536)
      dimension    amp(mmax),per(mmax)
!---------- 1. rescale input variables:
      a=x(1)*scal(1)
      b=x(2)*scal(2)
      nyp=2.*(t(2)-t(1))
      do j=1,nfou
         amp(j)=x(3*j)*scal(3*j)
!----------- scaling not changed generlization of scal necessary
         per(j)=x(3*j+1)*(50.-nyp)+nyp
      enddo
!---------- 2. compute chi**2
      sum=0.
      do i=1,ndata
         sum2=0.
         do j=1,nfou
            sum2=sum2+amp(j)*sin(2.*pi*(t(i)/per(j)+x(3*j+2)))
         enddo
         sum=sum+ ( (a*t(i)+b+sum2-f(i))/sigma )**2
      enddo
!---------- 3. define fitness
      fit1b=1./sum

      return
      end

