!****k* src/kpoints_3d.f90
!
! NAME
!  kpoints_3d
! 
! AUTHOR
!  W. Hergert
!
! HISTORY
!  * 28.03.2022 :old program checked
!
! USAGE 
!  Calculation of k-mesh in the irreducible part of the 3D-BZ
!  for bcc, fcc, and sc lattice, is used mainly for DOS 
!  calculation by means of the tetrahedron method (Dresden code)
!
!  q3bcc - kpoints bbc
!  q3fcc - kpoints fcc
!  q3sc  - kpoints bb
!  gltb  - tetrahedrons bcc
!  gltf  - tetrahedrons fcc
!  qkub  - tetrahedrons sc
!
! DESCRIPTION
!  These are the subprograms for the tetrahedron method, but only for
!  cubic BZ's. The programs are mainly an implementation of the former 
!  electronic structure group of TU Dresden (Ziesche, Lehmann, Wonn)
!
! INPUTS
!  1) kpoint generation
!   bcc:
!
!      NKP = (N+1)*(N+2)*(2*N+3)/6
!
!       N   |   1   2   3   4   5    6    7    8    9   10   11
!     --------------------------------------------------------
!       NKP |   5  14  30  55  91  140  204  285  385  505  650
!
!
!   fcc:
!
!     NUMBER OF K-POINTS :  NKP= 9*N**2+1+2*N*(8*N**2+7)/3
!
!      N  |     1   2    3    4    5    6    7
!      ---------------------------------------
!      NKP|    20  89  240  505  916 1505 2304
!
!
!   n     - subdivision of GX line
!   nkp   - maximum number of k-points
!   rgx   - length of GX line
!
!
!  2) tetrahedrons 
!
!    bcc:
!    N    : 2*N intervals on GAMMA-H
!    MI   : .GE. (5*N**3+3*N**2-2*N)/3
!
!    fcc
!    N    : 4*N INTERVALS ON LINE GAMMA-X
!    NT   : (80*N**3+12*N**2-2*N)/3
!
!   OUTPUT
!
!  1) kpoint generation
!                                                                                                                                                                                       
!   q     - array contains k-points
!   we    - array contains weights (weights are one!)   
!
!  2) tetrahedrons
!     
!      I(j) : index of the j-th corner of the tetrahedron
!
!
!  LITERATURE
!
!   Lehmann, G. and Taut, M., On the Numerical Calculation of the 
!   Density of States and Related Properties, physica status solidi (b),
!   54, 469-477 (1972).   
!
! SOURCE
!
MODULE kpoints_3D
PUBLIC :: q3bcc,q3fcc,q3sc

CONTAINS
!
! q3bcc
!
      SUBROUTINE q3bcc(n,nkp,rgx,q,we)
      USE typedef
      IMPLICIT NONE
      INTEGER(I4B)               :: n,nkp,n1,m,kz,ky,kx, &
                                    kxe,kye,kze
      REAL(DP)                   :: rgx,qx,qy,qz
      REAL(DP), DIMENSION(3,nkp) :: q
      REAL(DP), DIMENSION(nkp)   :: we
!
      INTENT (IN)    n,nkp,rgx
      INTENT (OUT)   q,we
!
      N1=2*N
      M=1
      KZE=N+1
      DO KZ=1,KZE
         KYE=N-KZ+2
         DO KY=1,KYE
            KXE=2*(N-KZ-KY)+5
            DO KX=1,KXE
               QX=(KX+KY+KZ-3.)/N1*RGX
               QY=(KY+KZ-2.)/N1*RGX
               QZ=(KZ-1.)/N1*RGX
               Q(1,M)=QX
               Q(2,M)=QY
               Q(3,M)=QZ
               WE(M)=1.
               M=M+1
            END DO
         END DO
      END DO
      END SUBROUTINE q3bcc
!
!  q3fcc
!
!
      SUBROUTINE q3fcc(n,nkp,rgx,q,we)
      USE typedef
      IMPLICIT NONE
      INTEGER(I4B)                 :: n,nkp,n1,m,kze,kye,kxe, &
                                      kx,ky,kz,n2
      REAL(DP)                     :: rgx,qx,qy,qz
      REAL(DP), DIMENSION(3,nkp)   :: q
      REAL(DP), DIMENSION(nkp)     :: we
!
      INTENT (IN)        n,nkp,rgx
      INTENT (OUT)       q,we
!
      N1=4*N
      M=1
      KZE=2*N+1
      DO KZ=1,KZE
         KYE=3*N+2-KZ-KZ/2
         DO KY=1,KYE
            KXE=N1-KZ-KY+3
            N2=6*N+6-3*KZ-2*KY
            IF(N2.LT.KXE)KXE=N2
            DO KX=1,KXE
               QY=(KX+KY+KZ-3.)/N1*RGX
               QX=(KY+KZ-2.)/N1*RGX
               QZ=(KZ-1.)/N1*RGX
               Q(2,M)=QY
               Q(1,M)=QX
               Q(3,M)=QZ
               WE(M)=1.D0
               M=M+1
           END DO
        END DO
      END DO
      END SUBROUTINE q3fcc
!
! q3sc
!
      SUBROUTINE q3sc(n,nst,rx,q,we)
      USE typedef
      IMPLICIT NONE
      INTEGER(I4B)                 :: n,nst,ns,i,l,k
      REAL(DP), DIMENSION(3,nst)   :: q
      REAL(DP), DIMENSION(nst)     :: we
      REAL(DP)                     :: rx,dx
!
      INTENT (IN)   n,nst,rx
      INTENT (OUT)  q,we
!
      DX=RX/(N-1)
      NS=1
      DO I=1,N
         DO L=I,N
            DO K=1,I
               Q(1,NS)=(I-1)*DX
               Q(2,NS)=(L-1)*DX
               Q(3,NS)=(K-1)*DX
               WE(NS)=1.D0
               NS=NS+1
            END DO
         END DO
      END DO
      END SUBROUTINE q3sc
!      
!  GLTB : tetrahedrons for bcc structure
!
      SUBROUTINE GLTB(I,MA,MI,N,NS)
      USE typedef
      IMPLICIT NONE
      INTEGER(I4B) :: m,n1,n2,n3,m1,m2,m3,m4,nt,n,l,ma,ns, &
                      k,nj,j,mi
      INTEGER(I4B), DIMENSION(4) :: i
      INTENT(IN)                 :: ma,mi,n,ns
      INTENT(OUT)                :: i
!
      m=mi
      M=1
      N1=1
      DO 100 L=1,N
         N2=N1+2*(N-L)+3
         N3=N1+(N-L+2)**2
         M1=N1
         M2=N2
         M4=N3
         N1=N3
         NT=2*(N-L)+1
         IF(NT.LT.1) GOTO 50
         DO 2 K=1,NT
            IF(M.NE.MA) GOTO 2
            I(1)=M1+K-1
            I(2)=M1+K
            I(3)=M2+K-1
            I(4)=M4+K-1
           RETURN
    2   M=M+1
50      IF(NT.LT.2) GOTO 60
        DO 3 K=2,NT
           IF(M.NE.MA) GOTO 10
           I(1)=M1+K-1
           I(2)=M2+K-1
           I(3)=M4+K-2
           I(4)=M4+K-1
           RETURN
10     M=M+1
       IF(M.NE.MA) GOTO 3
       I(1)=M1+K-1
       I(2)=M2+K-2
       I(3)=M2+K-1
       I(4)=M4+K-2
       RETURN
    3  M=M+1
60     IF(M.NE.MA) GOTO 11
       I(1)=M1+NT
       I(2)=M1+NT+1
       I(3)=M2+NT-1
       I(4)=M4+NT-1
       RETURN
11     M=M+1
       NJ=N-L
       IF(NJ.LT.1) GOTO 70
       DO 1 J=1,NJ
          M1=M2
          M3=M4
          M2=M2+NT
          M4=M4+NT
          NT=NT-2
          IF(NT.LT.1) GOTO 80
          DO 5 K=1,NT
             IF(M.NE.MA) GOTO 12
             I(1)=M1+K-1
             I(2)=M1+K
             I(3)=M2+K-1
             I(4)=M3+K
             RETURN
12        M=M+1
          IF(M.NE.MA) GOTO 5
          I(1)=M1+K
          I(2)=M3+K
          I(3)=M3+K+1
          I(4)=M4+K-1
          RETURN
    5     M=M+1
80        IF(M.NE.MA) GOTO 13
          I(1)=M1
          I(2)=M2
          I(3)=M3
          I(4)=M3+1
          RETURN
13        M=M+1
          IF(M.NE.MA) GOTO 14
          I(1)=M2
          I(2)=M3
          I(3)=M3+1
          I(4)=M4
          RETURN
14        M=M+1
          IF(M.NE.MA) GOTO 15
          I(1)=M1+NT
          I(2)=M2+NT-1
          I(3)=M3+NT+1
          I(4)=M4+NT-1
          RETURN
15        M=M+1
          IF(M.NE.MA) GOTO 16
          I(1)=M1+NT
          I(2)=M1+NT+1
          I(3)=M2+NT-1
          I(4)=M3+NT+1
          RETURN
16        M=M+1
          IF(NT.LT.2) GOTO 70
          DO 1 K=2,NT
             IF(M.NE.MA) GOTO 17
             I(1)=M1+K-1
             I(2)=M2+K-2
             I(3)=M2+K-1
             I(4)=M4+K-2
             RETURN
17           M=M+1
             IF(M.NE.MA) GOTO 18
             I(1)=M2+K-1
             I(2)=M3+K
             I(3)=M4+K-2
             I(4)=M4+K-1
             RETURN
18           M=M+1
             IF(M.NE.MA) GOTO 1
             I(1)=M1+K-1
             I(2)=M2+K-1
             I(3)=M3+K
             I(4)=M4+K-2
             RETURN
    1        M=M+1
70    CONTINUE
100   CONTINUE
      M=M-1
      END SUBROUTINE GLTB      
!
! GLTF : tetrahedrons for the fcc structure
!
      SUBROUTINE GLTF(I,MA,M,N,NS)
      USE typedef
      IMPLICIT NONE
      INTEGER(I4B) :: I(4),M1(20),M2(20),M3(20),MT(20)
      INTEGER(I4B) :: nx,ns,m,mc,n,mb,mh,ma,nk,j,l, &
                      n1,n2,n3,n4,n5,n6,n7,n8,k,k1,na,nt, &
		      nu1,nu2,nu,no
!
      nx=ns
      M=1
      MC=N
                   MB=N
      M1(1)=0
                M2(1)=4*N+1
                              M3(1)=(2*N+1)*(3*N+1)+N*N
                                                          MT(1)=4*N
      IF(N.LT.2) GOTO 786
      DO 1 L=2,N
      MH=MB
      MB=MC+1
      MC=MH
      M1(L)=M3(L-1)
                      M2(L)=M1(L)+MT(L-1)
      M3(L)=M3(L-1)+(2*N+3-2*L)*(3*N+1)+MC*MB
1     MT(L)=MT(L-1)-1
786   CONTINUE
      NU=N+1
               NO=2*N
      IF(MC.EQ.MB) GOTO 3
      MC=(3*N+2)/2
      MB=MC
                                                    GO TO 4
3     MC=(3*N+3)/2
      MB=MC-1
       IF(NO.LT.NU) GOTO 250
    4 DO 2 L=NU,NO
      M1(L)=M3(L-1)
                      M2(L)=M1(L)+MT(L-1)
      M3(L)=M3(L-1)+MC*MB
      MH=MB-1
      MB=MC-2
      MC=MH
    2 MT(L)=MT(L-1)-3
250    CONTINUE
      M3(2*N+1)=M3(2*N)+1
                            NU=3*N
                                     NU1=2
                                             NU2=1
      IF(NO.LT.1) GOTO 350
      DO 9 K=1,NO
      N1=M1(K)
                 N2=M2(K)
                            N3=M3(K)
                                       NT=MT(K)
      IF(NT.LT.1) GOTO 450
      DO 6 L=1,NT
      IF(M.NE.MA) GOTO 6
      I(1)=N1+L
      I(2)=N1+L+1
      I(3)=N2+L
      I(4)=N3+L
      RETURN
    6 M=M+1
450   IF(NT.LT.2) GOTO 550
      DO 7 L=2,NT
      IF(M.NE.MA) GOTO 200
      I(1)=N1+L
      I(2)=N2+L-1
      I(3)=N2+L
      I(4)=N3+L-1
      RETURN
200   M=M+1
      IF(M.NE.MA) GOTO 7
      I(1)=N1+L
      I(2)=N2+L
      I(3)=N3+L-1
      I(4)=N3+L
      RETURN
7     M=M+1
550   CONTINUE
              IF(K.LE.N) GO TO 15
      N5=N2-2
                N6=N2+NT
                           N7=N3+NT
      IF(M.NE.MA) GOTO 300
      I(1)=N5
      I(2)=N6
      I(3)=N6+1
      I(4)=N7
      RETURN
300   M=M+1
      IF(M.NE.MA) GOTO 40
      I(1)=N5
      I(2)=N5+1
      I(3)=N6+1
      I(4)=N7
      RETURN
40    M=M+1
      IF(M.NE.MA) GOTO 50
      I(1)=N5+1
      I(2)=N2
      I(3)=N6+1
      I(4)=N7
      RETURN
50    M=M+1
      IF(K.EQ.NO) GO TO 8
   15 NA=2*N-2*K+1
                     NU=NU-NU1
      MH=NU1
               NU1=NU2
                         NU2=MH
                                  N4=N3
      IF(NU.LT.1) GOTO 350
      DO 9 K1=1,NU
      N1=N2
              N2=N2+NT
                         N3=N4
                                 N4=N4+NT
                                            NT=NT-1
      IF(K1-NA)10,20,30
   30 N2=N2+1
   20 NT=NT-1
10    IF(NT.LT.1) GOTO 650
       DO 11 L=1,NT
      IF(M.NE.MA) GOTO 11
      I(1)=N1+L
      I(2)=N1+L+1
      I(3)=N2+L
      I(4)=N3+L+1
      RETURN
11    M=M+1
650   CONTINUE
      IF(M.NE.MA) GOTO 60
      I(1)=N1+1
      I(2)=N2+1
      I(3)=N3+1
      I(4)=N3+2
      RETURN
60    M=M+1
      IF(M.NE.MA) GOTO 70
      I(1)=N2+1
      I(2)=N3+1
      I(3)=N3+2
      I(4)=N4+1
      RETURN
70    M=M+1
      IF(NT.LT.2) GOTO 750
      DO 12 L=2,NT
      IF(M.NE.MA) GOTO 80
      I(1)=N1+L
      I(2)=N2+L-1
      I(3)=N2+L
      I(4)=N4+L-1
      RETURN
80    M=M+1
      IF(M.NE.MA) GOTO 90
      I(1)=N1+L
      I(2)=N3+L
      I(3)=N3+L+1
      I(4)=N4+L-1
      RETURN
90    M=M+1
      IF(M.NE.MA) GOTO 100
      I(1)=N2+L
      I(2)=N3+L+1
      I(3)=N4+L-1
      I(4)=N4+L
      RETURN
100   M=M+1
      IF(M.NE.MA) GOTO 12
      I(1)=N1+L
      I(2)=N2+L
      I(3)=N3+L+1
      I(4)=N4+L-1
      RETURN
   12 M=M+1
750   CONTINUE
              IF(K1.LT.NA) GO TO 9
      N5=N1+NT+1
                   N6=N2+NT
                              N7=N3+NT+1
                                           N8=N4+NT
      IF(K1.EQ.NA) GO TO 13
      IF(M.NE.MA) GOTO 110
      I(1)=N5+1
      I(2)=N5+2
      I(3)=N6+1
      I(4)=N7+1
      RETURN
110   M=M+1
13    IF(M.NE.MA) GOTO 120
      I(1)=N5
      I(2)=N5+1
      I(3)=N6+1
      I(4)=N7+1
      RETURN
120   M=M+1
      IF(M.NE.MA) GOTO 130
      I(1)=N5
      I(2)=N6
      I(3)=N6+1
      I(4)=N8
      RETURN
130   M=M+1
      IF(M.NE.MA) GOTO 140
      I(1)=N5
      I(2)=N7
      I(3)=N7+1
      I(4)=N8
      RETURN
140   M=M+1
      IF(M.NE.MA) GOTO 150
      I(1)=N5
      I(2)=N6+1
      I(3)=N7+1
      I(4)=N8
      RETURN
150   M=M+1
    9 CONTINUE
350   CONTINUE
    8 N3=M3(1)
      DO 14 K=1,N
      N1=N3
              N2=M3(2*K)-1
                             N3=M3(2*K+1)
      IF(M.NE.MA) GOTO 160
      I(1)=N1-3
      I(2)=N1
      I(3)=N2
      I(4)=N2+1
      RETURN
160   M=M+1
      IF(M.NE.MA) GOTO 170
      I(1)=N1-3
      I(2)=N1-2
      I(3)=N1
      I(4)=N2+1
      RETURN
170   M=M+1
      IF(M.NE.MA) GOTO 180
      I(1)=N1-2
      I(2)=N1-1
      I(3)=N1
      I(4)=N2+1
      RETURN
180   M=M+1
      IF(M.NE.MA) GOTO 14
      I(1)=N1
      I(2)=N2
      I(3)=N2+1
      I(4)=N3
      RETURN
   14 M=M+1
      M=M-1
END SUBROUTINE gltf
!      
! TETKUB : tetrahedrons sc lattice
!      
      SUBROUTINE TETKUB(ITET,NT,N)
      USE typedef
      IMPLICIT NONE 
      INTEGER(I4B) :: i,j,k,n1,n2,n3,n,nt,l,l1,l2, &
                      j1,j2,i1
      INTEGER(I4B) :: ITET(4,NT)
      J=0
      I=1
      K=1
   10 N1=J*(J+1)*(3*N-2*J+2)/6+1
      N3=(J+1)*(J+2)*(3*N-2*J)/6-J-1
      IF(I.EQ.N1)GOTO 1
! ZERLEGUNG DER PRISMEN MIT QUADRATISCHER GRUNDFLAECHE IN TETRAEDER
   12 IF(J.EQ.0)GOTO 5
      ITET(1,K)=I
      ITET(2,K)=(J+2)*(I-N1-J-1)/(J+1)+N3+J+2
      ITET(3,K)=ITET(2,K)+J+2
      ITET(4,K)=ITET(3,K)+1
      ITET(1,K+1)=I
      ITET(2,K+1)=I+1
      ITET(3,K+1)=ITET(2,K)+1
      ITET(4,K+1)=ITET(4,K)
      ITET(1,K+2)=I
      ITET(2,K+2)=ITET(2,K)
      ITET(3,K+2)=ITET(2,K)+1
      ITET(4,K+2)=ITET(4,K)
      ITET(1,K+3)=I
      ITET(2,K+3)=I+J+1
      ITET(3,K+3)=ITET(3,K)
      ITET(4,K+3)=ITET(4,K)
      ITET(1,K+4)=I
      ITET(2,K+4)=I+1
      ITET(3,K+4)=I+J+2
      ITET(4,K+4)=ITET(4,K)
      ITET(1,K+5)=I
      ITET(2,K+5)=I+J+1
      ITET(3,K+5)=I+J+2
      ITET(4,K+5)=ITET(4,K)
      IF(J.LE.1) GOTO 51
      I1=J-1
      DO 8 L=1,I1
      DO 7 L2=1,6
      DO 7 L1=1,4
      J1=K+L2-1
      J2=K+L2+5
      ITET(L1,J2)=ITET(L1,J1)+1
    7 CONTINUE
    8 K=K+6
   51 K=K+6
    5 ITET(1,K)=I+J
      ITET(2,K)=(J+2)*(I-N1-J-1)/(J+1)+N3+2*J+2
      ITET(3,K)=ITET(2,K)+J+2
      ITET(4,K)=ITET(3,K)-1-J
      ITET(1,K+1)=I+J
      ITET(2,K+1)=I+2*J+1
      ITET(3,K+1)=ITET(3,K)
      ITET(4,K+1)=ITET(3,K)+1
      ITET(1,K+2)=I+J
      ITET(2,K+2)=ITET(3,K)
      ITET(3,K+2)=ITET(4,K)
      ITET(4,K+2)=ITET(4,K+1)
      K=K+3
      I=I+J+1
      IF(I.LE.N3) GOTO 12
      J=J+1
      I=J+N3+1
      GOTO 10
    1 IF(J.EQ.0) GOTO 4
!   ZERLEGUNG DER PRISMEN MIT DREIECKIGER GRUNDFLAECHE IN TETRAEDER
      ITET(1,K)=I
      ITET(2,K)=I+1
      ITET(3,K)=I+J+2
      ITET(4,K)=N3+J+3
      ITET(1,K+1)=ITET(1,K)
      ITET(2,K+1)=ITET(3,K)-1
      ITET(3,K+1)=ITET(3,K)
      ITET(4,K+1)=ITET(4,K)-1
      ITET(1,K+2)=ITET(1,K)
      ITET(2,K+2)=ITET(3,K)
      ITET(3,K+2)=ITET(4,K)-1
      ITET(4,K+2)=ITET(4,K)
      IF(J.LE.1) GOTO 53
      I1=J-1
      DO 9 L=1,I1
      DO 3 L2=1,3
      DO 3 L1=1,4
      J1=K+L2-1
      J2=K+L2+2
      ITET(L1,J2)=ITET(L1,J1)+1
    3 CONTINUE
    9 K=K+3
   53 K=K+3
    4 ITET(1,K)=N1+J
      ITET(2,K)=N1+2*J+1
      ITET(3,K)=N3+2*J+2
      ITET(4,K)=ITET(3,K)+1
      IF(K.EQ.NT)GOTO 99
      K=K+1
      I=I+J+1
      GOTO 10
   99 CONTINUE
      RETURN
END SUBROUTINE tetkub
!
END MODULE kpoints_3d      
