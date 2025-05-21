!******* Function to round a real number to the nearest integer *****

      INTEGER FUNCTION IRND(A)
      IMPLICIT REAL*8 (A-H,O-Z)

      IF(A.GT.0) THEN
       K = INT(A)
       B = A - DFLOAT(K)
       IF(B.GT.0.5d0) K = K + 1
      ELSE
       K = INT(A)
       B = A - DFLOAT(K)
       IF(B.LT.-0.5d0) K = K - 1
      END IF
      IRND = K
      RETURN
      END
	  

C******* DIAGONALIZE THE REAL SYMMETRIC MATRIX *******************
C
C            ..... good for small matrices .....
C
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
      IMPLICIT REAL*8 (A-H,O-Z) 
      PARAMETER (NMAX=250)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.
11      CONTINUE
        V(IP,IP)=1.
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
	
      PAUSE '50 iterations should never happen'
      RETURN
      END
C
C**** SORT THE EIGENVALUES AND VECTORS OF THE DIAGONALIZED MATRIX *********
C
      SUBROUTINE EIGSRT(D,V,N,NP)
C
C  Modified by R. Krems to sort the eigenvalues and 
C  eigenvectors in the ascending rather then descending 
C  order. 
C               April 20, 2000, Goteborg
C 
C  Uses a selection sort algorithm.

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(NP),V(NP,NP)
      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).LE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END
C
C***** SORT THE ARRAY IN THE ASCENDING ORDER *************************
C
      SUBROUTINE SHELL(N,NDIM,ARR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ALN2I=1.4426950, TINY=1.E-5)
      DIMENSION ARR(NDIM)
      LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
      M=N
      DO 12 NN=1,LOGNB2
        M=M/2
        K=N-M
        DO 11 J=1,K
          I=J
3         CONTINUE
          L=I+M
          IF(ARR(L).LT.ARR(I)) THEN
            T=ARR(I)
            ARR(I)=ARR(L)
            ARR(L)=T
            I=I-M
            IF(I.GE.1)GO TO 3
          ENDIF
11      CONTINUE
12    CONTINUE
      RETURN
      END
C
C
C***** Computes the weights W(j) and Abscissas X(j) ******************
C***** for the Gauss-Hermite Quadrature Integration ******************
C
      SUBROUTINE gauher(x,w,n)
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(N),W(N)
      PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=100)
      m=(n+1)/2
      do 13 i=1,m
        if(i.eq.1)then
          z=sqrt(dfloat(2*n+1))-1.85575d0*(2*n+1)**(-.16667d0)
        else if(i.eq.2)then
          z=z-1.14d0*n**.426d0/z
        else if (i.eq.3)then
          z=1.86d0*z-.86d0*x(1)
        else if (i.eq.4)then
          z=1.91d0*z-.91d0*x(2)
        else
          z=2.d0*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=dsqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp
          if(dabs(z-z1).le.EPS) goto 1
12      continue
        pause 'too many iterations in gauher'
1       x(i)=z
        x(n+1-i)=-z
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
      return
      END
C
C*******  Abscissas for integration ************************************************
C
      SUBROUTINE ABSCISS(RMIN, RMAX, NUMP, XX)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION XX(NUMP)     
      RSTEP = (RMAX - RMIN)/DFLOAT(NUMP-1)
      RR = RMIN
      DO IP = 1, NUMP 
       XX(IP) = RR
       RR = RR + RSTEP
      END DO
      RETURN
      END
C
C***********  ASSOCIATED LEGENDRE POLYNOMIALS **********
C
C
C Computes the associated Legendre polynomial P^m_l(x)
C Here m and l are integers satisfying m=[0,l], x=[-1,1]
C
      REAL*8 FUNCTION plgndr(l,m,x)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      IF (m.lt.0.or.m.gt.l.or.DABS(x).gt.1.) THEN
      STOP 'bad arguments to plgndr'
      END IF
      pmm=1.D0   
C        ... Compute P^m_m
      IF(m.gt.0) THEN
      somx2=DSQRT((1.D0-x)*(1.D0+x))
      fact=1.D0
      DO i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2
      ENDDO
      ENDIF
C      
      IF(l.eq.m) THEN
          plgndr=pmm
      ELSE
        pmmp1=x*(2*m+1)*pmm   
C                 .....Compute P^m_{m+1}
       IF(l.eq.m+1) THEN
         plgndr=pmmp1
       ELSE                    
C           .... Compute P^m_l, l>m+1
        DO ll=m+2,l
           pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
           pmm=pmmp1
           pmmp1=pll
        ENDDO
        plgndr=pll
        ENDIF
       ENDIF
       RETURN 
       END 
C
C
C********* Spherical Harmonic *************************
C
       SUBROUTINE YJM(THETA,J,M,HARM)
C
C       ROUTINE TO EVALUATE THE SPHERICAL HARMONICS
C          Yjm(Theta,0)
C
C       INPUT:
C            Angle Theta (in rad), J and M values
C       OUTPUT:
C            The value of the spherical harmonic - HARM
C       USED FUNCTIONS:
C            Faq - to compute the factorial values
C            Plgndr - to compute associated Legendre polynomials
C
C                By Roman Krems,
C              June 15, 2000, Goteborg, 12.20 p.m.
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      SIGN = 1.d0
      IF(M.LT.0) THEN
       M = -M
       SIGN = -1.d0
       IF(MOD(M,2).EQ.0) SIGN = 1.d0
      END IF
C
      PI = DACOS(-1.d0)
      V1 = DFLOAT(J-M)
      V2 = DFLOAT(J+M)
C
      FACT = DFLOAT(2*J+1)/4.d0/PI*FAQ(V1)/FAQ(V2)
      FACT2 = DSQRT(FACT)
C
      COSANG = DCOS(THETA)       
C
      PL = PLGNDR(J,M,COSANG)
C
      HARM = SIGN*PL*FACT2
C
      RETURN
      END
C
C***************** Compute the factorial of v ******************
C
C
      REAL*8 FUNCTION faq(v)
C
C     *** calculates the factorial of z
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*8(I)
      PARAMETER (a0=1.0d0, a1=-0.5748646d0, a2=0.9512363,
     ¤  a3=-0.6998588, a4=0.4245549, a5=0.1010678)
C
       z = v
      IF (Z.GE.0.) then
       r = z - int(z)
       p = ((((a5*r+a4)*r + a3)*r +a2)*r + a1)*r + a0
       r = 1.d0
        s = z
        do while(s.gt.z-int(z)) 
          r = r*s
          s = s - 1
        end do
         q = r*p
      else
           if (z.gt.-1.) then
            r = (z+1) - int(z+1)
            p = ((((a5*r+a4)*r+a3)*r+a2)*r+a1)*r + a0
            q = p/(z+1)
           else
            q = 9.99999d9
           end if
      end if 
         faq = q
      RETURN
      END 
C
C
C
C****  Routine to generate the Abscissa and Weights  ******
C***       Gauss-Legendre Quadrature Integration   ********
C
C
      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (EPS=3.D-14)
      DIMENSION X(N), W(N)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
C
C
C
C************** 3J-SYMBOLS **********************************
C

      FUNCTION THREEJ (J1,J2,J3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     COMPUTATION OF SPECIAL WIGNER 3J COEFFICIENT WITH
C     VANISHING PROJECTIONS.  SEE EDMONDS, P. 50.
C
C     THIS VERSION EVALUATES BINOM AND PARITY IN-LINE
C       SHOULD IMPROVE EFFICIENCY, ESPECIALLY ON CRAY;
C       ALSO GIVES IMPROVEMENT ON AMDAHL  (SG: 20 DEC 92)
C
C     STATEMENT FUNCTION FOR DELTA ASSOCIATED W/ RACAH AND SIXJ SYMBOLS
C     DELTA(I,J,K)= DSQRT(1.D0/ ( BINOM(I+J+K+1,I+J-K) *
C    1                 BINOM(K+K+1,I-J+K) * DBLE(K+J-I+1) )  )
C
      I1=J1+J2+J3
      IF (I1-2*(I1/2).NE.0)    GO TO 8
    1 I2=J1-J2+J3
      IF (I2) 8,2,2
    2 I3=J1+J2-J3
      IF (I3) 8,3,3
    3 I4=-J1+J2+J3
      IF (I4) 8,4,4
    4 I5=I1/2
      I6=I2/2
      SIGN=1.D0
      IF (I5-2*(I5/2).NE.0) SIGN=-SIGN
C   7 THREEJ=SIGN*DELTA(J1,J2,J3)*BINOM(I5,J1)*BINOM(J1,I6)
C     B1,B2 ARE BINOM ASSOCIATED W/ DELTA
      N=J1+J2+J3+1
      M=J1+J2-J3
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B1=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 101 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  101   B = C/F
        B1 = B
      ENDIF
      N=J3+J3+1
      M=J1-J2+J3
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B2=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 102 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  102   B = C/F
        B2 = B
      ENDIF
      DELTA=DSQRT(1.D0/(B1*B2*(J3+J2-J1+1)))
C     B3=BINOM(I5,J1),  B4=BINOM(J1,I6)
      N=I5
      M=J1
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B3=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 103 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  103   B = C/F
        B3 = B
      ENDIF
      N=J1
      M=I6
      NM = N-M
      MNM = MIN(NM,M)
      IF(MNM.LE.0) THEN
        B4=1.D0
      ELSE
        FN = N+1
        F = 0.D0
        B = 1.D0
        DO 104 I = 1,MNM
        F = F+1.D0
        C = (FN-F)*B
  104   B = C/F
        B4 = B
      ENDIF
      THREEJ=SIGN*DELTA*B3*B4
      RETURN
    8 THREEJ=0.D0
      RETURN
      END
C
C**********************************************************************
C
C
C**********************************************************************
C
C
C--------------------------------------------------------------------
      FUNCTION XF6J(A,B,E,D,C,F)
*
*                                    | A  B  E |
*   PROGRAM TO COMPUTE THE 6J SYMBOL {         }
*                                    | D  C  F |
*   AUTHOR: B. FOLLMEG
*   CURRENT REVISION DATE: 4-MAY-1997
*
* -------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FCT/ SI(100000)
      DATA TOL,ZERO,ONE,TWO /1.D-10,0.D0,1.D0,2.D0/
      X=ZERO
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( A B E)
      IF ((E .GT. (A + B)) .OR. (E .LT. ABS(A - B))) GOTO 40
      SUM = A + B + E
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IABE = NINT(SUM)
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( D C E)
      IF ((E .GT. (C + D)) .OR. (E .LT. ABS(C - D))) GOTO 40
      SUM = D + C + E
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IDCE = NINT(SUM)
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( A C F)
      IF ((F .GT. (A + C)) .OR. (F .LT. ABS(A - C))) GOTO 40
      SUM = A + C + F
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IACF = NINT(SUM)
* CHECK TRIANGULAR CONDITIONS FOR TRIAD ( D B F)
      IF ((F .GT. (D + B)) .OR. (F .LT. ABS(D - B))) GOTO 40
      SUM = D + B + F
      IF (MOD(SUM,ONE) .GT. TOL) GOTO 40
      IDBF = NINT(SUM)
      IABDC = NINT(A + B + D + C)
      IAEDF = NINT(A + E + D + F)
      IBECF = NINT(B + E + C + F)
      MINCHI = MAX(IABE,IDCE,IACF,IDBF,-1)
      MAXCHI = MIN(IABDC,IAEDF,IBECF) - MINCHI
* INDICES FOR DELTAS
      DELTA = ZERO
      I2A = NINT(TWO*A) - 1
      I2B = NINT(TWO*B) - 1
      I2C = NINT(TWO*C) - 1
      I2D = NINT(TWO*D) - 1
      I2E = NINT(TWO*E) - 1
      I2F = NINT(TWO*F) - 1
* DELTA(ABE)
      J1 = IABE - I2A
      J2 = IABE - I2B
      J3 = IABE - I2E
      J4 = IABE + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
* DELTA(DCE)
      J1 = IDCE - I2D
      J2 = IDCE - I2C
      J3 = IDCE - I2E
      J4 = IDCE + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
* DELTA(ACF)
      J1 = IACF - I2A
      J2 = IACF - I2C
      J3 = IACF - I2F
      J4 = IACF + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
* DELTA(DBF)
      J1 = IDBF - I2D
      J2 = IDBF - I2B
      J3 = IDBF - I2F
      J4 = IDBF + 2
      DELTA = DELTA + SI(J1) + SI(J2) + SI(J3) - SI(J4)
      DELTA = 0.5D0 * DELTA
      IABDC = IABDC - MINCHI
      IAEDF = IAEDF - MINCHI
      IBECF = IBECF - MINCHI
      IABE = MINCHI - IABE
      IDCE = MINCHI - IDCE
      IACF = MINCHI - IACF
      IDBF = MINCHI - IDBF
      ABDC = DBLE(IABDC - MAXCHI)
      AEDF = DBLE(IAEDF - MAXCHI)
      BECF = DBLE(IBECF - MAXCHI)
      ABE = DBLE(MAXCHI + IABE + 1)
      DCE = DBLE(MAXCHI + IDCE + 1)
      ACF = DBLE(MAXCHI + IACF + 1)
      DBF = DBLE(MAXCHI + IDBF + 1)
* LOOP OVER CHI
      X = 1.D0
      IPOWER = 0
      IF (MAXCHI .LE. 0) GOTO 30
      II = MINCHI + MAXCHI + 2
      DO 10 ICHI = 1, MAXCHI
         XCHI = DBLE(ICHI)
         XA = (ABDC + XCHI) * (AEDF + XCHI) * (BECF + XCHI)
         XB = (ABE - XCHI) * (DCE - XCHI) * (ACF - XCHI) * (DBF - XCHI)
         X = 1.D0 - XA * (II - ICHI) * X / XB
10    CONTINUE
      IF (X) 20,40,30
20    X = -X
      IPOWER = 1
30    X = DLOG(X) + SI(MINCHI+2)-SI(IABDC+1)-SI(IAEDF+1)-SI(IBECF+1)
     :            - SI(IABE+1)-SI(IDCE+1)-SI(IACF+1)-SI(IDBF+1) + DELTA
      IPOWER = IPOWER + MINCHI
      X = (-1.D0)**IPOWER * DEXP(X)
40    XF6J = X
      RETURN
      END
C
      SUBROUTINE FACTORIAL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/FCT/FACT(100000)
      FACT(1)=0.D0
      FACT(2)=0.D0
      N=10000
 2    DO I=3,N
         FACT(I)=FACT(I-1)+DLOG(DFLOAT(I-1))
      ENDDO
      RETURN
      END
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C

C**********************************************************************
C
C
C
      FUNCTION THRJ(F1,F2,F3,G1,G2,G3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     SMALL CHANGES 31 JUL 95 (SG)
      SAVE MUNG,X,Y
      PARAMETER (MXIX=302)
      DIMENSION X(MXIX),Y(MXIX)
      DATA MUNG/0/
      IF (MUNG.EQ.21) GO TO 69
      MUNG = 21
      X(1) = 0.D0
      DO 100 I = 1, MXIX-1
      A = I
      X(I+1) = LOG(A) +X(I)
      Y(I+1) = LOG(A)
  100 CONTINUE
   69 IF(F1-ABS(G1)) 1,13,13
   13 IF(F2-ABS(G2))1,14,14
   14 IF(F3-ABS(G3))1,15,15
   15 SUM=F1+F2+F3
      NSUM=SUM+.001D0
      IF(SUM-NSUM)2,2,1
    1 THRJ=0.D0
      RETURN
    2 IF(ABS(G1+G2+G3)-1.D-08)3,3,1
    3 IF(F1+F2-F3)1,4,4
    4 IF(F1+F3-F2)1,5,5
    5 IF(F2+F3-F1)1,6,6
    6 J1=2.D0*F3+2.001D0
      J2=F1+F2-F3+1.001D0
      J3=F1-F2+F3+1.001D0
      J4=-F1+F2+F3+1.001D0
      J5=F1+F2+F3+2.001D0
      J6=F1+G1+1.001D0
      J7=F1-G1+1.001D0
      J8=F2+G2+1.001D0
      J9=F2-G2+1.001D0
      J10=F3+G3+1.001D0
      J11=F3-G3+1.001D0
      IF(J5.GT.MXIX) THEN
        WRITE(6,601) J5,MXIX
  601   FORMAT(' *** DIMENSION ERROR IN THRJ - INDEX.GT.MXIX',2I5)
        STOP
      ENDIF
      R=0.5D0*(Y(J1)+X(J2)+X(J3)+X(J4)-X(J5)
     1+X(J6)+X(J7)+X(J8)+X(J9)+X(J10)+X(J11))
      SUM=0.D0
      F=-1
      KZ=-1
    7 KZ=KZ+1
      F=-F
      J1=KZ+1
      J2=F1+F2-F3-KZ+1.001D0
      IF(J2)20,20,8
    8 J3=F1-G1-KZ+1.001D0
      IF(J3)20,20,9
    9 J4=F2+G2-KZ+1.001D0
      IF(J4)20,20,10
   10 J5=F3-F2+G1+KZ+1.001D0
      IF(J5)7,7,11
   11 J6=F3-F1-G2+KZ+1.001D0
      IF(J6)7,7,12
   12 JMAX=MAX(J1,J2,J3,J4,J5,J6)
      IF(JMAX.GT.MXIX) THEN
        WRITE(6,601) JMAX,MXIX
        STOP
      ENDIF
      S=-(X(J1)+X(J2)+X(J3)+X(J4)+X(J5)+X(J6))
      SUM=SUM+F*EXP(R+S)
      GO TO 7
   20 INT=ABS(F1-F2-G3)+0.0001D0
      VAL=((-1.D0)**INT)*SUM/SQRT(2.D0*F3+1.D0)
      IF(ABS(VAL).LE.1.D-6) VAL=0.D0
      THRJ=VAL
      RETURN
      END
C

C
C***********************************************************************
C*************  Invert a matrix A(ND,ND) ******************
C
C
      SUBROUTINE INV (A,N,ND)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ND,ND)
      DIMENSION PIVOT(1000000),IPIVOT(1000000),INDEX(1000000,2)
      IF(N.EQ.1) GO TO 741
   15 DO 20 J=1,N
   20 IPIVOT(J) = 0
   30 DO 550 I= 1,N
   40 AMAX =0.0
   45 DO 105 J=1,N
   50 IF(IPIVOT(J)-1) 60,105,60
   60 DO 100 K=1,N
   70 IF(IPIVOT(K)-1) 80,100,740
   80 IF(DABS(AMAX)-DABS(A(J,K))) 85,100,100
   85 IROW=J
   90 ICOLUM=K
   95 AMAX=A(J,K)
  100 CONTINUE
  105 CONTINUE
  110 IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
  130 IF(IROW-ICOLUM) 150,260,150
  150 DO 200 L=1,N
  160 SWAP= A(IROW,L)
  170 A(IROW,L) = A(ICOLUM,L)
  200 A(ICOLUM,L) = SWAP
  260 INDEX (I,1) =IROW
  270 INDEX(I,2)=ICOLUM
  310 PIVOT(I) =A(ICOLUM,ICOLUM)
  330 A(ICOLUM,ICOLUM) =1.00D0
  340 DO 350  L=1,N
  350  A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
  380 DO 550 L1 =1,N
  390 IF(L1-ICOLUM) 400,550,400
  400 T =A(L1,ICOLUM)
  420 A(L1,ICOLUM) = 0.0D0
  430 DO 450  L=1,N
  450 A(L1,L) = A(L1,L)-A(ICOLUM,L)*T
  550 CONTINUE
  600 DO 710 I=1,N
  610 L = N+1-I
  620 IF(INDEX(L,1)-INDEX(L,2)) 630,710,630
  630 JROW = INDEX (L,1)
  640 JCOLUM = INDEX (L,2)
  650 DO  705 K= 1,N
  660 SWAP = A(K,JROW)
  670 A(K,JROW) = A(K,JCOLUM)
  700 A(K,JCOLUM) = SWAP
  705 CONTINUE
  710 CONTINUE
  740 RETURN
  741 A(1,1) = 1.D0/A(1,1)
      RETURN
       END
C
C**************  MULTIPLY MATRICES ***********************
C
      SUBROUTINE MULT1 (A,B,C,NV,NO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NV,NV),B(NV,NV),C(NV,NV)
      DO 1 I1=1,NO
      DO 1 I2=1,NO
      R1=0.D0
      DO 2 J=1,NO
2     R1=R1+B(I1,J)*C(J,I2)
1     A(I1,I2)=R1
      RETURN
      END
C
C**************  MULTIPLY MATRICES ***********************
C
      SUBROUTINE MULTTR (A,B,C,NV,NO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NV,NV),B(NV,NV),C(NV,NV)
      DO 1 I1=1,NO
      DO 1 I2=1,NO
      R1=0.D0
      DO 2 J=1,NO
2     R1=R1+B(I1,J)*C(I2,J)
1     A(I1,I2)=R1
      RETURN
      END
C
C**************  MULTIPLY MATRICES ***********************
C
      SUBROUTINE TRMULT (A,B,C,NV,NO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NV,NV),B(NV,NV),C(NV,NV)
      DO 1 I1=1,NO
      DO 1 I2=1,NO
      R1=0.D0
      DO 2 J=1,NO
2     R1=R1+B(J,I1)*C(J,I2)
1     A(I1,I2)=R1
      RETURN
      END
C
C*************************** MULTS ***********************
C
      SUBROUTINE MULTS (A,B,C,ND,NV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ND,ND), C(ND,ND)
      DO 1 I1=1,NV
      DO 1 I2=1,NV
 1    A(I1,I2)=C(I1,I2)*B
      RETURN
      END
C
C************************** ADDS *************************
C
      SUBROUTINE ADDS (A,B,C,ND,NV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ND,ND), C(ND,ND)
      DO 1 I=1,NV
    1 A(I,I)=C(I,I)+B
      RETURN
      END
C
C************************** MULTT ************************
C
      SUBROUTINE MULTT (NV,L,M,N,A,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NV,NV),B(NV,NV),C(NV,NV)
      DO 1 I1=1,L
      DO 1 I2=1,N
      R1=0.D0
      DO 2 J=1,M
2     R1=R1+B(I1,J)*C(J,I2)
1     A(I1,I2)=R1
      RETURN
      END
C
C************************** ADD ***************************
C
      SUBROUTINE ADD (A,B,C,ND,NV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(ND,ND),B(ND,ND),C(ND,ND)
      DO 1 I1=1,NV
      DO 1 I2=1,NV
    1 A(I1,I2)=B(I1,I2)+C(I1,I2)
      RETURN
      END
C
C ************* BESSEL FUNCTIONS ***************************
C
      SUBROUTINE RBES(N,Z,ZJN,ZYN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA ONE/1.0D0/
      DATA TWO /2.0D0/
      DATA HALF /0.5D0/
      DATA ZERO /0.0D0/
      DATA THREE /3.0D0/
      SINZ=DSIN(Z)
      COSZ=DCOS(Z)
      ZINV=ONE/Z
      AJ=SINZ
      BJ=SINZ*ZINV-COSZ
      AY=-COSZ
      BY=-SINZ-COSZ*ZINV
      IF(N-1) 1, 2,3
1     ZJN=AJ
      ZYN=AY
      RETURN
2     ZJN=BJ
      ZYN=BY
      RETURN
3     XN=N
      DELF=ZINV+ZINV
      FACTOR=DELF+ZINV
      IF(Z.GT.XN) GO TO 100
      ZSQ=Z**2
      TOP=XN+XN
      XJ=THREE
      C=ZSQ/THREE
25    ZYN=BY*FACTOR-AY
      AY=BY
      BY=ZYN
      FACTOR=FACTOR+DELF
      XJ=XJ+TWO
      C=C*(Z/XJ)
      IF(XJ.LT.TOP) GO TO 25
      FACTOR=-HALF*ZSQ
      U=ONE
      TERM=ONE
      DU=ONE
      DTERM=ONE
      D2NP3=TOP+THREE
      DEN=D2NP3
      DFACT=ZERO
35    DFACT=DFACT+ONE
      TERM=TERM*(FACTOR/(DFACT*DEN))
      U=U+TERM
      DEN=DEN+TWO
      DTERM=DTERM*(FACTOR/(DFACT*DEN))
      DU=DU+DTERM
      IF(DABS(TERM).GT.1.0D-15) GO TO 35
      ZJN=U*C
      RETURN
100   CONST=ZINV*XN
      TOP=CONST+CONST
200   AJ=FACTOR*BJ-AJ
      AY=FACTOR*BY-AY
      FACTOR=FACTOR+DELF
      IF(FACTOR.GT.TOP) GO TO 250
      BJ=FACTOR*AJ-BJ
      BY=FACTOR*AY-BY
      FACTOR=FACTOR+DELF
      IF(FACTOR.LT.TOP) GO TO 200
      ZJN=BJ
      ZYN=BY
      RETURN
250   ZJN=AJ
      ZYN=AY
      RETURN
      END
C
C************************************************************  
C      Spherical Bessel and Neumann functions,
C     i.e. Bessel functions of non-integer order
C
C************************************************************
C
C
      SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)
      IMPLICIT REAL*8 (A-H,O-Z)
CU    USES bessjy
      PARAMETER (RTPIO2=1.2533141d0)
      if(n.lt.0.or.x.le.0.) then
        print*, ' wrong values', n, x 
        pause 'bad arguments in sphbes'
      end if
      order=n+0.5
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.*x)
      syp=factor*ryp-sy/(2.*x)
      return
      END                  
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (EPS=1.e-10,FPMIN=1.e-30,MAXIT=10000,XMIN=2.,
     *PI=3.141592653589793d0)
CU    USES beschb
      if(x.le.0..or.xnu.lt.0.) pause 'bad arguments in bessjy'
      if(x.lt.XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        h=del*h
        if(d.lt.0.d0)isign=-isign
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessjy; try asymptotic expansion'
1     continue
      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do 12 l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
12    continue
      if(rjl.eq.0.d0)rjl=EPS
      f=rjpl/rjl
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2).lt.EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del).lt.(1.d0+abs(sum))*EPS)goto 2
13      continue
        pause 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do 14 i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di).lt.FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci).lt.FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli).lt.EPS)goto 3
14      continue
        pause 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do 15 i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
15    continue
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
      END
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NUSE1=5,NUSE2=5)
CU    USES chebev
      DIMENSION c1(7), c2(8)
      DATA c1/-1.142022680371172d0, 6.516511267076d-3, 3.08709017308d-4,
     *-3.470626964d-6, 6.943764d-9, 3.6780d-11, -1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6, -3.3126120d-8, 2.42310d-10, -1.70d-13, -1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END      
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
      FUNCTION chebev(a,b,c,m,x)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(m)
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END
C
CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
       DOUBLE PRECISION FUNCTION 
     *         XF9J(F1, F2, F3, 
     *               F4, F5, F6, 
     *               F7, F8, F9)
       IMPLICIT NONE 
       DOUBLE PRECISION :: F1, F2, F3, F4, F5, F6, F7, F8, F9
       DOUBLE PRECISION :: X, Xmin, Xmax, Sum
       INTEGER :: NumX, I
       INTEGER :: IRND
       DOUBLE PRECISION :: XF6J
C 
       Xmax = 0.d0
       IF((F1+F9).GT.Xmax) Xmax = F1 + F9
       IF((F4+F8).GT.Xmax) Xmax = F4 + F8
       IF((F2+F6).GT.Xmax) Xmax = F2 + F6
       Xmin = 1000.d0
       IF(DABS(F1-F9).LT.Xmin) Xmin = DABS(F1 - F9)
       IF(DABS(F4-F8).LT.Xmin) Xmin = DABS(F4 - F8)
       IF(DABS(F2-F6).LT.Xmin) Xmin = DABS(F2 - F6)
C
       NumX = IRND(Xmax - Xmin) + 1
       Sum = 0.d0
C
        X = Xmin - 1.d0
        DO I = 1, NumX 
        X = X + 1.d0
        Sum = Sum + 
     *      (-1.d0)**IRND(2.d0*X)*
     *       (2.d0*X+1.d0)*
     *     XF6J(F1,F4,F7,F8,F9,X)*
     *     XF6J(F2,F5,F8,F4,X,F6)*
     *     XF6J(F3,F6,F9,X,F1,F2)
        END DO
       XF9J = Sum
       RETURN
       END FUNCTION XF9J
C
C**********************************************************************
C
      DOUBLE PRECISION FUNCTION CG(F1,F2,F3,
     *                             G1,G2,G3)
      IMPLICIT NONE
      DOUBLE PRECISION :: F1,F2,F3, G1,G2,G3
      INTEGER :: IRND
      DOUBLE PRECISION :: THRJ
C
      CG = (-1.d0)**IRND(F1-F2+G3)*DSQRT(2.d0*F3+1.d0)*
     *  THRJ(F1,F2,F3,G1,G2,-G3)
      RETURN
      END 
C
C**********************************************************************

	SUBROUTINE SELECTION_SORT(NumberOfElements,Array)
!	Sort an double precision array in ascending order, using a selection algorithm	

	INTEGER :: I, NumberOfElements, LocationSmallest
	INTEGER :: MINLOC_array(1)
	DOUBLE PRECISION :: SmallestElement, Array(NumberOfElements)
	
	
	DO I=1,NumberOfElements-1
	  
	  SmallestElement = MINVAL(Array(I:NumberOfElements))
	  MINLOC_array = (MINLOC(Array(I:NumberOfElements)))
	  LocationSmallest = (I - 1) + MINLOC_array(1)
	  Array(LocationSmallest) = Array(I)
	  Array(I) = SmallestElement
	  	 
	END DO
	
	END 




