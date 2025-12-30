!****m* src/glbi.f90
!
! NAME
!  glbi.f90
! 
! AUTHOR
!  Theory group of Dresden University
!
! HISTORY
!  * 28.03.2022 checked
!
! USAGE 
!  Calculation of density of states by means of tetrahedron method
!  
!
! DESCRIPTION
!  This is a very old version  version of the tetrahedron method. The 
!  method
!  was developed by the the theory group in Dresden.
!
! INPUTS
!   EMIN - minimmum of energy axis
!     DE - energy step width
!   EMAX - maximum of energy axis
!     EF -
!  SPLIT - if energies in a tetrahedron are equal, they will 
!          be split to avoid division by zero
!  FNORM - normlization factor
!   NDIM - number of energy points
!    NST - number of k-points
!   NPGA - subdivision of irred. Part of BZ 
!   NPGM - subdivision of irred. part of BZ 
!    DEL - ?
!      A - array of weight factors for PDOS calculation for a single band
!      E - array  of energies in one band
!      Q - array of k-vectors
!    TET - subroutine to describe tetrahedrons EXTERNAL)
!
! OUTPUT
!    D - density of states  
!
! WARNINGS
!   q2dsq stops if number of k-points to be calculated is larger than nkp
!
! SOURCE
!
      SUBROUTINE GLBI(EMIN,DE,EMAX,EF,SPLIT,FNORM,NDIM,                     &
                      NST,NPGM,NPGA,DEL,A,E,Q,D,TET)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /TRANS/ DELT(4),AT(4),QT(3,4),EMINC,DEC,EMAXC,SPLITC
    DIMENSION ET(4),DEL(NST),A(NST),E(NDIM),Q(3,NST),ITET(4),D(NDIM),     &
                   DELT1(4),AT1(4),QT1(3,4)               
      EXTERNAL TET
      DATA IP/0/
!
      EMINC =EMIN
      DEC   =DE
      EMAXC =EMAX
      SPLITC=SPLIT
      NENIN =(EMAX-EMIN)/DE+1.1
!
      IF(NDIM.GT.NENIN) GOTO 222
         WRITE(IP,200)
200      FORMAT('0WRONG DIMENSIONS')
         RETURN
!
222   DO 1 I=1,1000000
         CALL TET(ITET,I,NT,NPGM,NPGA)
         DO 2 K=1,4
            IF(I.GT.NT) RETURN
            I1=ITET(K)
            IF(I1.GT.NST)RETURN
            AT(K)=A(I1)
            ET(K)=E(I1)
            DELT(K)=DEL(I1)
            DO 21 L=1,3
               QT(L,K)=Q(L,I1)
21          CONTINUE
2        CONTINUE 
!
! ASSIGNMENT OF VALUES FOR ONE TETRAHEDRON
!
         DO 3 K=1,3
            LO=4-K
            DO 31 L=1,LO
               IF(ET(L+1).GE.ET(L))GOTO 3
               AT1(1)=AT(L)
               AT(L)=AT(L+1)
               AT(L+1)=AT1(1)
               ET1=ET(L)
               ET(L)=ET(L+1)
               ET(L+1)=ET1
               DELT1(1)=DELT(L)
               DELT(L)=DELT(L+1)
               DELT(L+1)=DELT1(1)
               DO 32 M=1,3
                  QT1(1,1)=QT(M,L)
                  QT(M,L)=QT(M,L+1)
                  QT(M,L+1)=QT1(1,1)
32             CONTINUE
31          CONTINUE
3        CONTINUE
!
! ORDERING OF ET(I) COUPLED WITH AT(I),DELT(I) AND QT(I,K)
!
         DO 41 K=1,3
            IF((ET(K+1)-ET(K)).LT.SPLIT) ET(K+1)=ET(K)+SPLIT
41       CONTINUE
!
! SPLITTING OF ET(I)
!
         IF(EF.LE.ET(1)) GOTO 1
!
! EMPTY TETRAHEDRON
!
         IF(EF-ET(4)) 6,5,5
5        CALL GLOT(D,FNORM,NDIM)
         GOTO 1
!
! FULL TETRAHEDRON
!
6        IF(EF-ET(2)) 7,7,8
7        DO 71 K=2,4
            X=(EF-ET(1))/(ET(K)-ET(1))
            DELT(K)=DELT(1)+X*(DELT(K)-DELT(1))
            AT(K)=AT(1)+X*(AT(K)-AT(1))
            DO 72 L=1,3
               QT(L,K)=QT(L,1)+X*(QT(L,K)-QT(L,1))
72          CONTINUE
71       CONTINUE
         CALL GLOT(D,FNORM,NDIM)
         GOTO 1
!
! FULL EDGE OF THE TETRAHEDRON
!
8        IF(EF-ET(3)) 10,9,9
9        CALL GLOT(D,FNORM,NDIM)
         DO 91 K=1,3
            X=(EF-ET(4))/(ET(K)-ET(4))
            DELT(K)=DELT(4)+X*(DELT(K)-DELT(4))
            AT(K)=-(AT(4)+X*(AT(K)-AT(4)))
            DO 92 L=1,3
               QT(L,K)=QT(L,4)+X*(QT(L,K)-QT(L,4))
92          CONTINUE
91       CONTINUE
         AT(4)=-AT(4)
         CALL GLOT(D,FNORM,NDIM)
         GOTO 1
!
! EMPTY EDGE OF THE TETRAHEDRON
!
10       DO 101 K=3,4
            X=(EF-ET(1))/(ET(K)-ET(1))
            DELT1(K)=DELT(1)+X*(DELT(K)-DELT(1))
            AT1(K)=AT(1)+X*(AT(K)-AT(1))
            DO 102 L=1,3
               QT1(L,K)=QT(L,1)+X*(QT(L,K)-QT(L,1))
102         CONTINUE
101      CONTINUE
!
         DO 103 K=3,4
            X=(EF-ET(2))/(ET(K)-ET(2))
            DELT(K)=DELT(2)+X*(DELT(K)-DELT(2))
            AT(K)=AT(2)+X*(AT(K)-AT(2))
            DO 104 L=1,3
               QT(L,K)=QT(L,2)+X*(QT(L,K)-QT(L,2))
104         CONTINUE
103      CONTINUE
!
         CALL GLOT(D,FNORM,NDIM)
         DELT(2)=DELT1(3)
         AT(2)=AT1(3)
         DO 105 L=1,3
            QT(L,2)=QT1(L,3)
105      CONTINUE
         CALL GLOT(D,FNORM,NDIM)
!
         DELT(3)=DELT1(4)
         AT(3)=AT1(4)
         DO 106 L=1,3
            QT(L,3)=QT1(L,4)
106      CONTINUE
         CALL GLOT(D,FNORM,NDIM)
1     CONTINUE
      WRITE(IP,202)
202   FORMAT('0FEHLER IN ''TET''')
      RETURN
      END
!
      SUBROUTINE GLOT(D,FNORM,NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION D(NDIM)
      COMMON /TRANS/ DELT(4),AT(4),QT(3,4),EMIN,DE,EMAX,SPLIT
      COMMON /INTER/ DEL(4),A(4),Q(9)
!
      TRI(PRO,SUM,DELX,AX,E,VOL)=VOL/2.0/PRO*(AX*(E-DELX)**2+SUM/3.0*(E-DELX)**3)
!
      NENIN=(EMAX-EMIN)/DE+1.1
      DO 9 I=1,4
         DEL(I)=DELT(I)
         A(I)=AT(I)
9     CONTINUE
!
      DO 11 K=1,3
         IO=4-K
         DO 10 I=1,IO
            IF(DEL(I+1).GE.DEL(I)) GOTO 10
            B=DEL(I)
            DEL(I)=DEL(I+1)
            DEL(I+1)=B
            C=A(I)
            A(I)=A(I+1)
            A(I+1)=C
10       CONTINUE
11    CONTINUE
!
! ORDERING OF DEL(I) COUPLED WITH A(I)
!
      DO 110 I=1,3
         DO 111 K=1,3
            L=3*(I-1)+K
            Q(L)=QT(K,I)-QT(K,4)
111      CONTINUE
110   CONTINUE
      VOL=Q(1)*Q(5)*Q(9) &
         +Q(4)*Q(8)*Q(3) &
         +Q(7)*Q(2)*Q(6) &
         -Q(3)*Q(5)*Q(7) &
         -Q(6)*Q(8)*Q(1) &
         -Q(9)*Q(2)*Q(4)
      VOL=ABS(VOL)
!
! VOLUME OF THE TETRAHEDRON
!
      DO 1 I=1,3
         IF((DEL(I+1)-DEL(I)).LT.SPLIT) DEL(I+1)=DEL(I)+SPLIT
1     CONTINUE
!
! SPLITTING OF DEL(I)
!
      NU=(DEL(1)-EMIN)/DE+2.000001
      NO=(DEL(2)-EMIN)/DE+1.000001
      IF(NU.GT.NO)    GOTO 3
      IF(NU.GT.NENIN) GOTO 7
      CALL SUMPRO(1,SUM1,PRO1)
      DO 2 K=NU,NO
         IF(K.GT.NENIN) GOTO 3
         IF(K.LT.1)     GOTO 2
         E=EMIN+(K-1)*DE
         D(K)=D(K)+TRI(PRO1,SUM1,DEL(1),A(1),E,VOL)*FNORM
2     CONTINUE
!
3     NU=(DEL(2)-EMIN)/DE+2.000001
      NO=(DEL(3)-EMIN)/DE+1.000001
      IF(NU.GT.NO) GOTO 5
      CALL SUMPRO(1,SUM1,PRO1)
      CALL SUMPRO(2,SUM2,PRO2)
      DO 4 K=NU,NO
         IF(K.GT.NENIN) GOTO 5
         IF(K.LT.1)     GOTO 4
         E=EMIN+(K-1)*DE
         D(K)=D(K)+TRI(PRO1,SUM1,DEL(1),A(1),E,VOL)*FNORM    &
                  -TRI(PRO2,SUM2,DEL(2),A(2),E,VOL)*FNORM
4     CONTINUE
!
5     NU=(DEL(3)-EMIN)/DE+2.000001
      NO=(DEL(4)-EMIN)/DE+1.000001
      IF(NU.GT.NO) GOTO 7
      CALL SUMPRO(4,SUM4,PRO4)
      DO 6 K=NU,NO
         IF(K.GT.NENIN) GOTO 7
         IF(K.LT.1)     GOTO 6
         E=EMIN+(K-1)*DE
         D(K)=D(K)+TRI(PRO4,SUM4,DEL(4),A(4),E,VOL)*FNORM
6     CONTINUE
7     CONTINUE
!
      RETURN
      END
!
      SUBROUTINE SUMPRO(L,SUM,PR0)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /INTER/ DEL(4),A(4),Q(9)
      PR0=1.0
      SUM=0.0
      DO 1 I=1,4
         IF(I.EQ.L) GO TO 1
         SUM=SUM+(A(I)-A(L))/(DEL(I)-DEL(L))
         PR0=PR0*(DEL(I)-DEL(L))
    1 CONTINUE
      IF(L.NE.1) PR0=-PR0
      RETURN
      END
