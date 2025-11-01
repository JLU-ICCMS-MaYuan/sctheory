      program xpk1b
c==================================================================
c     Driver program for non-linear least-squares problem
c     (Sect. 5.3)  
c     modified for input 
c==================================================================
      implicit  none
      common/data/ f(200),t(200),sigma,ndata,m
      integer   n, seed, i, status,ndata,m
      parameter (n=17)
      real      ctrl(12), x(n), fit, fit1b,x1(n),nyp
      real      f,t,sigma
      external  fit1b
c
c     First, initialize the random-number generator
c
      seed=13579
      call rninit(seed)
c
c     Initializations
c
      call finit
c
c     Set control variables (use defaults except for population size)
c      
      do 10 i=1,12
         ctrl(i) = -1
   10 continue
      ctrl(1)=50

c     Now call pikaia
      call pikaia(fit1b,n,ctrl,x,fit,status)
c
c     Print the results
      write(*,*) ' status: ',status
      write(*,*) '      x: '
      WRITE(*,50) x(1),x(2)
50    FORMAT(3x,5(F12.6,2x))
      DO i=1,5
         WRITE(*,50) x(3*i),x(3*i+1),x(3*i+2)
      ENDDO           
      write(*,*) '      f: ',fit
      write(*,20) ctrl
 20   format(   '    ctrl: ',6f11.6/10x,6f11.6)
      nyp=2.*(t(2)-t(1))
      x1(1)=10.*x(1)
      x1(2)=100.*x(2)
      DO i=1,5
         x1(3*i)=x(3*i)*100.
         x1(3*i+1)=x(3*i+1)*(50.-nyp)+nyp
         x1(3*i+2)=x(3*i+2)
      ENDDO
      WRITE(*,50) x1(1),x1(2)
      DO i=1,5
         WRITE(*,50) x1(3*i),x1(3*i+1),x1(3*i+2)
      ENDDO     
      OPEN(14,FILE='fir.parm',FORM='formatted',STATUS='unknown')
      DO i=1,17
         WRITE(14,30) i,x1(i)
      ENDDO
30    FORMAT(I5,2x,F14.6)            
c
      end
c***************************************************************
      subroutine finit
c===============================================================
c     Reads in synthetic dataset (see Figure 5.4)
c===============================================================
      implicit     none
      common/data/ f(200),t(200),sigma,ndata,m
      dimension    f0(200),vdum(11)
      real         f,t,f0,vdum,sigma,delt
      integer      ndata,m,i

      open(1,file='test1n.dat',form='formatted')

      read(1,*) ndata
      DO i=1,ndata
         read(1,*) t(i),f(i)
      ENDDO
c     Use 5 Fourier modes for the fit
      m=5

c     same error bar for all point
      sigma=5.
      delt=t(3)-t(1)

      return
      end
      
