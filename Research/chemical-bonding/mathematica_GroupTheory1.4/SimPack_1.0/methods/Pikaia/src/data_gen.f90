PROGRAM data_gen
IMPLICIT NONE
INTEGER it,i,m,k,n
REAL f(200),t(200),a(3),p(3),fi(3)
REAL tx,fx,a0,b,pi,su
!
DATA a/12.,10.,8./
DATA p/20.,9.,7.5/
DATA fi/0.5,0.8,0.8/
!
WRITE(6,*) 'Which test *'
READ(5,*) it
IF((it.EQ.1).OR.(it.EQ.2)) THEN
m=3
pi=3.1516926536
a0=1.0
b=20.0
n=200
DO i=1,200
   tx=(i-1)*0.5
   t(i)=tx
   su=0.0 
   DO k=1,m
      su=su+a(k)*SIN(2*pi*(tx/p(k)+fi(k))) 
   ENDDO
   f(i)=a0*tx+b+su
ENDDO
WRITE(8,1) n
DO i=1,n
   WRITE(8,2) t(i),f(i)
ENDDO 
ENDIF   
1 FORMAT(I5)
2 FORMAT(E14.6,2x,E14.6)
END PROGRAM data_gen
