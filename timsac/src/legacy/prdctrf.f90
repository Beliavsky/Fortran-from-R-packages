! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine prdctrf(n,p,q,h,d,k,l,jsw,yy,b,a,ww,s,y,yori,yd,x,z1,&
&z2,z3,zz1,zz2,zz3)
  use timsac_kinds, only: dp
  implicit none
!
!c      PROGRAM PRDCTR
!     PROGRAM 74.3.1.  PREDICTION PROGRAM
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "TIMSAC-74 A TIME SERIES ANALYSIS AND CONTROL PROGRAM PACKAGE(2
!         BY H. AKAIKE, E. ARAHATA AND T. OZAKI, COMPUTER SCIENCE MONOGRA
!         NO.6 MARCH 1976, THE INSTITUTE OF STATISTICAL MATHEMATICS
!-----------------------------------------------------------------------
!     THIS PROGRAM OPERATES ON A REAL RECORD OF A VECTOR PROCESS
!     Y(I) (I=1,N) AND COMPUTES PREDICTED VALUES. ONE STEP AHEAD
!     PREDICTION STARTS AT TIME P AND ENDS AT TIME Q. PREDICTION IS
!     CONTINUED WITHOUT NEW OBSERVATIONS UNTIL TIME Q+H.
!     BASIC MODEL IS THE AUTOREGRESSIVE MOVING AVERAGE
!     MODEL OF Y(I) WHICH IS GIVEN BY
!     Y(I)+B(1)Y(I-1)+...+B(K)Y(I-K) = X(I)+A(1)X(I-1)+...+A(L)X(I-L).
!
!    THE FOLLOWING INPUTS ARE REQUIRED:
!     (N,P,Q,H):
!                  N, LENGTH OF DATA
!                  P, ONE STEP AHEAD PREDICTION STARTING POSITION
!                  Q, LONG RANGE FORECAST STARTING POSITION
!                  H, MAXIMUM SPAN OF FORECAST (LESS THAN OR EQUAL TO 100)
!                  (Q+H MUST BE LESS THAN 1001)
!     JSW: JSW=0 FOR DIRECT LOADING OF AR-MA COEFFICIENTS,
!                  THE OUTPUTS OF PROGRAM MARKOV WITH ICONT=0.
!            JSW=1 FOR LOADING OF THE OUTPUTS OF PROGRAM MARKOV,
!                  THE OUTPUTS OF PROGRAM MARKOV WITH ICONT=1.
!     (D,K,L):
!              D, DIMENSION OF THE VECTOR Y(I)
!              K, AR-ORDER (LESS THAN OR EQUAL TO 10)
!              L, MA-ORDER (LESS THAN OR EQUAL TO 10)
!     N,L,K,H,P,Q,D,JSW,ARE ALL INTEGERS
!     (DFORM(I),I=1,20): INPUT FORMAT STATEMENT IN ONE CARD,
!                             FOR EXAMPLE, (8F10.4)
!     (NAME(I,J),I=1,20,J=1,D): NAME OF THE I-TH COMPONENT
!     (Y(I,J),I=1,N;J=1,D): ORIGINAL DATA
!     (B(I1,I2,J),I1=1,D,I2=1,D,J=1,K): AR-COEFFICIENT MATRICES.
!     FOR JSW=0,
!         (A(I1,I2,J),I1=1,D,I2=1,D,J=1,L): MA-COEFFICIENT MATRICES.
!     FOR JSW=1,
!         (W(I1,I2,J),I1=1,D,I2=1,D,J=1,L): IMPULSE RESPONSE MATRICES.
!     (S(I,J),I=1,D,J=1,D): INNOVATION VARIANCE MATRIX
!
!     THE OUTPUTS OF THIS PROGRAM ARE THE REAL
!     AND PREDICTED VALUES OF Y(I).
!
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      INTEGER H,H1,P,Q,D
!c      INTEGER DFORM(20),XX(121)
!c      INTEGER IKA/1HA/, IKB/1HB/, IKW/1HW/
!c      REAL*8 A,B,W,S,Z,C,CX
!c      REAL*8 CST0,CST1
!c      DIMENSION A(10,10,10),B(10,10,10),S(10,10),W(10,10,101)
!c      DIMENSION X(1000,10),C(10,10),SD(10,101),CX(10),EY(10),NAME(20,10)
!c      DIMENSION Y(1000,10),YORI(500,10)
!c      DIMENSION DMXT(10),DMIT(10),AV(10)
!xx      DIMENSION A(D,D,L),B(D,D,K),S(D,D),WW(D,D,L),W(D,D,H+1)
!xx      DIMENSION X(N,D),C(D,D),SD(D,H+1),CX(D),EY(D)
!xx      DIMENSION YY(N,D),YORI(H+1,D)
!xx      DIMENSION DMXT(D),DMIT(D),AV(D)
!c      EQUIVALENCE(W(1,1,1),X(1,1))
!xx      DIMENSION Y(Q+H,D),YD(Q+H,D)
!xx      DIMENSION Z1(Q+H,D),Z2(Q+H,D),Z3(Q+H,D)
!xx      DIMENSION ZZ1(Q+H,D),ZZ2(Q+H,D),ZZ3(Q+H,D)
integer n, p, q ,h, d, k, l, jsw
real(dp) yy(n,d), b(d,d,k), a(d,d,l), ww(d,d,l), s(d,d),&
&y(q+h,d), yori(h+1,d), yd(q+h,d), x(n,d),&
&z1(q+h,d), z2(q+h,d), z3(q+h,d), zz1(q+h,d),&
&zz2(q+h,d), zz3(q+h,d)
! local
integer i, i2, ii, ij, im1, iqh, isr, isw, ix, iy, j, jj, jx, jy,&
&jz, kk, h1
real(dp) w(d,d,h+1), c(d,d),sd(d,h+1), cx(d), ey(d),&
&dmxt(d), dmit(d), av(d), cst0, cst1, z, ccc, an,&
&ave, dmax, dmin, yyd
!
!x      INTEGER*1  TMP(1)
!x      CHARACTER  CNAME*80
!
!      DATA A/1000*0.0D-00/, B/1000*0.0D-00/, W/10000*0.0D-00/
!c      DATA A/1000*0.0D-00/, B/1000*0.0D-00/, W/10100*0.0D-00/
!c      DATA S/100*0.0D-00/, C/100*0.0D-00/
!c      DATA CX/10*0.0D-00/
!c      DATA X/10000*0.0/,SD/1010*0.0/,EY/10*0.0/,Y/10000*0.0/
!c      DATA YORI/5000*0.0/,DMXT/10*0.0/,DMIT/10*0.0/,AV/10*0.0/
!c      DATA X/10000*0.0D-00/,SD/1010*0.0D-00/,EY/10*0.0D-00/
!c      DATA Y/10000*0.0D-00/,YORI/5000*0.0D-00/
!c      DATA DMXT/10*0.0/,DMIT/10*0.0/,AV/10*0.0/
!c      INTEGER XX(121)
!c      DATA XX/121*1H /
!c      DATA K1 / 1H* /, K2 / 1HX /, K3 / 1H+ /, K4 / 1HY /, K5 / 1H  /
!c      DATA K6/1H!/
!xx      CHARACTER KSTOR
!
!     INPUT / OUTPUT DATA FILE OPEN
!c      CALL SETWND
!c      CALL FLOPN2(NFL)
!c      IF (NFL.EQ.0) GO TO 999
!x      IER=0
!x      LU=3
!x      DO 100 I = 1,80
!x  100 CNAME(I:I) = ' '
!x      I = 1
!x      IFG = 1
!x      DO WHILE( (IFG.EQ.1) .AND. (I.LE.80) )
!x         IF ( TMP(I).NE.ICHAR(' ') ) THEN
!x            CNAME(I:I) = CHAR(TMP(I))
!x            I = I+1
!x         ELSE
!x            IFG = 0
!x         END IF
!x      END DO
!x      IF ( I.GT.1 ) THEN
!x         IFG = 1
!x         OPEN (LU,FILE=CNAME,IOSTAT=IVAR)
!x         IF (IVAR .NE. 0) THEN
!xcx            WRITE(*,*) ' ***  prdctr temp FILE OPEN ERROR :',CNAME,IVAR
!x            IER=IVAR
!x            IFG=0
!x         END IF
!x      END IF
!
!c      MJ1=10
!c      MJ2=10
cst0=0.0d-00
cst1=1.0d-00
isw=1
h1=h+1
iqh=q+h
!
!      DATA INITIALIZE
!
if (jsw.ne.0) a(1:d,1:d,1:l) = cst0
if (jsw.eq.0) ww(1:d,1:d,1:l) = cst0
w(1:d,1:d,1:h1) = cst0
c(1:d,1:d) = cst0
cx(1:d) = cst0
x(1:n,1:d) = cst0
sd(1:d,1:h1) = cst0
ey(1:d) = cst0
y(1:iqh,1:d) = cst0
yori(1:h1,1:d) = cst0
dmxt(1:d) = cst0
dmit(1:d) = cst0
av(1:d) = cst0
!
!     INITIAL CONDITION INPUT
!c      READ(5,800) N,D,P,Q,H
!c      READ(5,800) JSW
!c      DO 199 J=1,D
!c      READ(5,801) (NAME(I,J),I=1,20)
!c  199 CONTINUE
!c      READ(5,801) (DFORM(I),I=1,20)
!     ORIGINAL DATA VECTOR (Y(I),I=1,N) INPUT
!c      DO 100 J=1,D
!c      READ(5,DFORM) (Y(I,J),I=1,N)
!c  100 CONTINUE
!xx      DO 70 I = 1,N
do 71 i = 1,n
do 70 j = 1,d
y(i,j) = yy(i,j)
70 continue
71 continue
!
!c      READ(5,800) D,K,L
!c      WRITE(6,900)
!c      WRITE(6,901) N,D,K,L,P,Q,H,JSW
!c      WRITE(6,902) (DFORM(I),I=1,20)
!c      DO 1101 J=1,D
!c      WRITE(6,903) J,(NAME(I,J),I=1,20)
!c      WRITE(6,904) (Y(I,J),I=1,N)
!c 1101 CONTINUE
!
!     AR-COEFFICIENTS INPUT
!c      CALL REMT3X(B,D,D,K,ISW,MJ1,MJ1,MJ2,IKB)
!c      IF  (JSW.NE.0) GO TO 80
!     MA-COEFFICIENTS INPUT
!c      CALL REMT3X(A,D,D,L,ISW,MJ1,MJ1,MJ2,IKA)
!c      GO TO 81
!     IMPULSE RESPONSE INPUT
!c   80 CALL REMT3X(W,D,D,L,ISW,MJ1,MJ1,MJ2,IKW)
do 80 ii = 1,d
!xx      DO 80 JJ = 1,D
!xx      DO 80 KK = 1,L
do 79 jj = 1,d
do 78 kk = 1,l
w(ii,jj,kk) = ww(ii,jj,kk)
78 continue
79 continue
80 continue
!     INNOVATION VARIANCE INPUT
!     COMMON SUBROUTINE CALL
!c   81 CALL REMATX(S,D,D,ISW,MJ1,MJ1)
!c      WRITE(6,1357)
!c      CALL SUBMPR(S,D,D,MJ1,MJ1)
if  (jsw.eq.0) go to 99
!
!     MA-COEFFICIENTS COMPUTATION
!     A(M)=W(M)+B(1)W(M-1)+...+B(M-1)W(1)+B(M)
!     M=1,L
do  84 i=1,l
do  83 jx=1,d
do  82 jy=1,d
!xx   82 A(JX,JY,I)=W(JX,JY,I)+B(JX,JY,I)
a(jx,jy,i)=w(jx,jy,i)+b(jx,jy,i)
82 continue
83 continue
84 continue
do  89 i=2,l
im1=i-1
do  88 ix=1,im1
iy=i-ix
do  87 jx=1,d
do  86 jy=1,d
z=cst0
do  85 jz=1,d
!xx   85 Z=Z+B(JX,JZ,IX)*W(JZ,JY,IY)
z=z+b(jx,jz,ix)*w(jz,jy,iy)
85 continue
!xx   86 A(JX,JY,I)=A(JX,JY,I)+Z
a(jx,jy,i)=a(jx,jy,i)+z
86 continue
87 continue
88 continue
89 continue
do 8210 i=1,l
!xx      DO 8220 IX=1,D
do 8221 ix=1,d
do 8220 jx=1,d
!xx 8220 C(IX,JX)=A(IX,JX,I)
c(ix,jx)=a(ix,jx,i)
8220 continue
8221 continue
!c      WRITE(6,8121) I
!c      CALL SUBMPR(C,D,D,MJ1,MJ1)
8210 continue
!
!     W(I) CLEAR
!xx      DO  98 I=1,L
!xx      DO  97 JX=1,D
!xx      DO  96 JY=1,D
!xx   96 W(JX,JY,I)=CST0
!xx   97 CONTINUE
!xx   98 CONTINUE
w(1:d,1:d,1:l)=cst0
!
!     IMPULSE RESPONSE W(M),M=0,H COMPUTATION
!     W(0) = UNIT MATRIX
!     W(M) = A(M)-(B(1)W(M-1)+...+B(K)W(M-K)) (FOR M LESS THAN OR EQUAL
!     W(M) =-B(1)W(M-1)-...-B(K)W(M-K) (FOR M GREATER THAN L)
99 do  101 i=1,d
!xx  101 W(I,I,1)=CST1
w(i,i,1)=cst1
101 continue
do 120 i=2,h1
do  115 j=1,k
ij=i-j
!c      IF  (IJ) 115,115,102
if  (ij.le.0) go to 115
!xx  102 DO  105 JX=1,D
do  105 jx=1,d
do  104 ix=1,d
z=cst0
do  103 i2=1,d
!xx  103 Z=Z+(B(IX,I2,J)*W(I2,JX,IJ))
z=z+(b(ix,i2,j)*w(i2,jx,ij))
103 continue
w(ix,jx,i) = w(ix,jx,i) - z
104 continue
105 continue
115 continue
!c      IF  (I-L-1) 116,116,120
if  (i-l-1.gt.0) go to 120
!xx  116 DO  118 JX=1,D
do  118 jx=1,d
do  117 ix=1,d
!xx  117 W(IX,JX,I) = W(IX,JX,I)+A(IX,JX,I-1)
w(ix,jx,i) = w(ix,jx,i)+a(ix,jx,i-1)
117 continue
118 continue
120 continue
!
!     PREDICTION ERROR VARIANCE C(I) AND STANDARD DEVIATION SD(IX,I)
!     COMPUTATION (I=1,H1)
!     C(I)=W(0)*S*W(0)'+...+W(I)*S*W(I)' AND SD(IX,I) IS THE VECTOR OF
!     THE POSITIVE SQUARE ROOTS OF THE DIAGONAL ELEMENTS OF THE PREDICTI
!     ERROR VARIANCE MATRIX C(I).  ONLY THE DIAGONAL ELEMENTS OF C(I)
!     ARE COMPUTED.
do 140 i=1,h1
do  125 ix=1,d
do  124 jx=1,d
z=cst0
do 123 i2=1,d
!xx  123 Z=Z+W(IX,I2,I)*S(I2,JX)
z=z+w(ix,i2,i)*s(i2,jx)
123 continue
c(ix,jx)=z
124 continue
125 continue
do  139 ix=1,d
z=cst0
do  129 i2=1,d
!xx  129 Z=Z+C(IX,I2)*W(IX,I2,I)
z=z+c(ix,i2)*w(ix,i2,i)
129 continue
cx(ix)=cx(ix)+z
ccc=cx(ix)
!c      IF  (CCC) 130,130,131
!c  130 SD(IX,I)=0.0
if  (ccc.gt.0) go to 131
!xx  130 SD(IX,I)=0.0D-00
sd(ix,i)=0.0d-00
go to 139
!c  131 SD(IX,I)=SQRT(CCC)
131 sd(ix,i)=dsqrt(ccc)
139 continue
140 continue
!
!     SUBTRACTION OF THE MEAN VALUES FROM THE ORIGINAL DATA Y(I)
an=n
!c	AN=1.0/AN
an=1.0d-00/an
do  146  i=1,d
z=cst0
do  141  j=1,n
!xx  141 Z=Z+Y(J,I)
z=z+y(j,i)
141 continue
ave=z*an
av(i)=ave
do  145  j=1,n
y(j,i)=y(j,i)-ave
145 continue
146 continue
!
!
!     PREDICTIONS AND INNOVATIONS (X(I)) COMPUTATION OF Y(I) (I=1,Q+H)
!     FOR I GREATER THAN OR EQUAL TO Q, X(I) IS SET EQUAL TO 0.
!
!
do 6300 j=1,d
dmxt(j)=y(1,j)
!xx 6300 DMIT(J)=Y(1,J)
dmit(j)=y(1,j)
6300 continue
isr=0
do 300 i=1,iqh
!     EY, PREDICTED VALUE OF Y, COMPUTATION
do  153 j=1,d
!c  153 EY(J)=0.0
!xx  153 EY(J)=0.0D-00
ey(j)=0.0d-00
153 continue
!     B(1)Y(I-1)+...+B(K)Y(I-K) COMPUTATION
do  160 j=1,k
ij=i-j
!c      IF  (IJ) 160,160,154
if  (ij.le.0) go to 160
!xx  154 DO  159 IX=1,D
do  159 ix=1,d
z=cst0
do  158 i2=1,d
!xx  158 Z=Z+B(IX,I2,J)*Y(IJ,I2)
z=z+b(ix,i2,j)*y(ij,i2)
158 continue
ey(ix)=ey(ix)-z
159 continue
160 continue
!     A(1)X(I-1)+...+A(L)X(I-L) COMPUTATION
do  170 j=1,l
ij=i-j
!c      IF  (IJ) 170,170,161
!c  161 IF  (IJ-Q) 162,170,170
if  (ij.le.0) go to 170
!xx  161 IF  (IJ-Q.GE.0) GO TO 170
!xx  162 DO  169 IX=1,D
if  (ij-q.ge.0) go to 170
do  169 ix=1,d
z=cst0
do  168 i2=1,d
!xx  168 Z=Z+(A(IX,I2,J)*X(IJ,I2))
z=z+(a(ix,i2,j)*x(ij,i2))
168 continue
ey(ix)=ey(ix)+z
169 continue
170 continue
!     MAXIMUM, MINIMUM SEARCH
do  249 j=1,d
dmax=dmxt(j)
dmin=dmit(j)
!c      YD=Y(I,J)
!c      CALL MAXMIN(DMAX,DMIN,YD)
!c      YD=EY(J)
!c      CALL MAXMIN(DMAX,DMIN,YD)
yyd=y(i,j)
call maxmin(dmax,dmin,yyd)
yyd=ey(j)
call maxmin(dmax,dmin,yyd)
dmxt(j)=dmax
dmit(j)=dmin
249 continue
!     INNOVATION X(I) IS GIVEN BY (Y-EY).
!     IF I IS GREATER THAN OR EQUAL TO Q, THEN X(I)=0.
!c      IF(I-Q) 171,200,200
if(i-q.ge.0) go to 200
!xx  171 DO  172  J=1,D
do  172  j=1,d
x(i,j)=y(i,j)-ey(j)
172 continue
go to 300
200 isr=isr+1
do 250 j=1,d
yori(isr,j)=y(i,j)+av(j)
y(i,j)=ey(j)
250 continue
300 continue
!
do 281 i=1,iqh
do  280 j=1,d
!xx  280 Y(I,J)=Y(I,J)+AV(J)
y(i,j)=y(i,j)+av(j)
280 continue
281 continue
!
!
!     ********************
!     PRINT OUT
!     ********************
i2=0
!c	CALL HEADPR(1,P,Q,D,N)
do 500 ii=1,iqh
i=ii
!c      IF  (I-Q)   301,449,450
!c  301 IF  (I-P)   302,399,400
!c      IF  (I-Q)	  301,450,450
!c  301 IF  (I-P)	  500,400,400
if  (i-q.ge.0) go to 450
!xx  301 IF  (I-P.LT.0) GO TO 500
if  (i-p.lt.0) go to 500
!c  302 DO  310 J=1,D
!c      IF  (J-1) 303,303,304
!c  303 WRITE(6,910) I,Y(I,J)
!c      GO TO 310
!c  304 WRITE(6,911) Y(I,J)
!c  310 CONTINUE
!c      WRITE(6,920)
!c      GO TO 500
!c  399 CALL HEADPR(I,P,Q,D,N)
!xx  400 DO  410 J=1,D
do  410 j=1,d
!c      YD=Y(I,J)-X(I,J)
yd(i,j)=y(i,j)-x(i,j)
!c      IF  (J-1) 403,403,404
!c  403 WRITE(6,910) I,Y(I,J),YD,X(I,J)
!c      GO TO 410
!c  404 WRITE(6,911) Y(I,J),YD,X(I,J)
410 continue
!c      WRITE(6,920)
go to 500
!c  449 CALL HEADPR(I,P,Q,D,N)
450 i2=i2+1
do  480 j=1,d
!c      YD=Y(I,J)
!c      Z1=YD+SD(J,I2)
!c      Z2=Z1+SD(J,I2)
!c      Z3=Z2+SD(J,I2)
yd(i,j)=y(i,j)
z1(i,j)=yd(i,j)+sd(j,i2)
z2(i,j)=z1(i,j)+sd(j,i2)
z3(i,j)=z2(i,j)+sd(j,i2)
!      IF  (J-1) 453,453,454
!  453 WRITE(6,912) I,YD,Z1,Z2,Z3
!      IF(I.LE.N) WRITE(6,937) YORI(I2,J)
!      GO TO 455
!  454 WRITE(6,913) YD,Z1,Z2,Z3
!      IF(I.LE.N) WRITE(6,937) YORI(I2,J)
!c      IF  ((J-1).LE.0) THEN
!c         IF (I.LE.N) THEN
!c            WRITE(6,9370) I,YORI(I2,J),YD,Z1,Z2,Z3
!c         ELSE
!c            WRITE(6,912) I,YD,Z1,Z2,Z3
!c         END IF
!c      ELSE
!c         IF (I.LE.N) THEN
!c            WRITE(6,9371) YORI(I2,J),YD,Z1,Z2,Z3
!c         ELSE
!c            WRITE(6,913) YD,Z1,Z2,Z3
!c         END IF
!c      END IF
!c  455 Z1=YD-SD(J,I2)
!c      Z2=Z1-SD(J,I2)
!c      Z3=Z2-SD(J,I2)
!xx  455 ZZ1(I,J)=YD(I,J)-SD(J,I2)
zz1(i,j)=yd(i,j)-sd(j,i2)
zz2(i,j)=zz1(i,j)-sd(j,i2)
zz3(i,j)=zz2(i,j)-sd(j,i2)
!c      WRITE(6,914) Z1,Z2,Z3
480 continue
!c      WRITE(6,920)
500 continue
!
!
!     ********************
!     GRAPHIC PRINT OUT
!     ********************
!x      IF (IFG .NE. 0) THEN
!x      KSTOR=K2
!x      DO  600  J=1,D
!x      DMAX=DMXT(J)
!x      DMIN=DMIT(J)
!c      FMAX=ABS(DMAX)
!c      FMIN=ABS(DMIN)
!c      FMAX=AMAX1(FMAX,FMIN)
!x      FMAX=DABS(DMAX)
!x      FMIN=DABS(DMIN)
!x      FMAX=DMAX1(FMAX,FMIN)
!x      FMIN=-FMAX
!xcc      YST=FMAX/60.0
!x      YST=FMAX/60.0D-00
!x      TFMIN=FMIN+AV(J)
!x      TFMID=AV(J)
!x      TFMAX=FMAX+AV(J)
!c      WRITE(6,931) J,(NAME(I,J),I=1,20)
!c      WRITE(6,932)
!c      WRITE(6,933)
!c      WRITE(6,934) TFMIN,TFMID,TFMAX
!c      WRITE(6,935)
!x      WRITE(LU,931) J
!x      WRITE(LU,932)
!x      WRITE(LU,933)
!x      WRITE(LU,934) TFMIN,TFMID,TFMAX
!x      WRITE(LU,935)
!
!x      DO 599 I=1,IQH
!x      YA=Y(I,J)-AV(J)
!x      XX(61)=K6
!c      IF  (I-Q)	 501,569,570
!c  501 IF  (I-P)	 502,550,550
!x      IF  (I-Q.EQ.0) GO TO 569
!x      IF  (I-Q.GT.0) GO TO 570
!x  501 IF  (I-P.GE.0) GO TO 550
!
!     REAL DATA
!x  502 CALL SBSCAL(YA,YST,IX)
!x      XX(IX)=K1
!c      WRITE(6,936) I,(XX(I2),I2=1,121)
!x      WRITE(LU,936) I,(XX(I2),I2=1,121)
!x      XX(IX)=K5
!x      GO TO 599
!
!     REAL DATA
!x  550 CALL SBSCAL(YA,YST,IX)
!     ONE-STEP PREDICTION
!x      YB=YA-X(I,J)
!x      CALL SBSCAL(YB,YST,JX)
!
!c  566 IF  (IX-JX) 568,567,568
!x  566 IF  (IX-JX.NE.0) GO TO 568
!x  567 XX(IX)=K3
!c      WRITE(6,936) I,(XX(I2),I2=1,121)
!x      WRITE(LU,936) I,(XX(I2),I2=1,121)
!x      XX(IX)=K5
!x      GO TO 599
!
!x  568 XX(IX)=K1
!x      XX(JX)=K2
!c      WRITE(6,936) I,(XX(I2),I2=1,121)
!x      WRITE(LU,936) I,(XX(I2),I2=1,121)
!x      XX(IX)=K5
!x      XX(JX)=K5
!x      GO TO 599
!
!x  569 KSTOR=K2
!x      K2=K4
!x      ISR=0
!
!     LONG RANGE PREDICTION
!x  570 CALL SBSCAL(YA,YST,JX)
!c  576 IF  (I-N) 577,577,598
!x  576 IF  (I-N.GT.0) GO TO 598
!     REAL DATA
!x  577 ISR=ISR+1
!x      YC=YORI(ISR,J)-AV(J)
!x      CALL SBSCAL(YC,YST,IX)
!x      GO TO 566
!
!x  598 XX(JX)=K4
!c      WRITE(6,936) I,(XX(I2),I2=1,121)
!x      WRITE(LU,936) I,(XX(I2),I2=1,121)
!x      XX(JX)=K5
!x  599 CONTINUE
!x      K2=KSTOR
!x  600 CONTINUE
!x      END IF
!
!c      CALL FLCLS2(NFL)
!c  999 CONTINUE
!x      IF (IFG.NE.0) CLOSE(LU)
!xx  800 FORMAT(8I5)
!xx  801 FORMAT(20A4)
!xx  900 FORMAT(1H ,' PROGRAM 74.3.1. PREDICTION')
!xx  901 FORMAT(/1H ,' INITIAL CONDITION: N=',I4,', D=',I2,', K=',I2,
!xx     A', L=',I4,', P=',I4,', Q=',I4,', H=',I4,', JSW=',I1)
!xx  902 FORMAT(/1H ,' ORIGINAL DATA (',20A4,'):')
!xx  903 FORMAT(1H ,' ** Y ',I2,2X,20A4)
!xx  904 FORMAT(1H ,10E12.5)
!c  910 FORMAT(1H ,'   N=',I5,4X,E20.5,4(4X,E20.5))
!c  911 FORMAT(1H ,14X,E20.5,4(4X,E20.5))
!c  912 FORMAT(1H ,'   N=',I5,24X,4(4X,E20.5))
!c  913 FORMAT(1H ,34X,4(4X,E20.5))
!c  914 FORMAT(1H ,58X,3(4X,E20.5))
!xx  910 FORMAT(1H ,'   N=',I5,4X,D20.5,4(4X,D20.5))
!xx  911 FORMAT(1H ,14X,D20.5,4(4X,D20.5))
!xx  912 FORMAT(1H ,'   N=',I5,24X,4(4X,D20.5))
!xx  913 FORMAT(1H ,34X,4(4X,D20.5))
!xx  914 FORMAT(1H ,58X,3(4X,D20.5))
!xx  920 FORMAT(1H ,4X)
!xx  931 FORMAT(/1H ,10X,'J=',I2,2X,20A4)
!xx  932 FORMAT(/1H ,50X,'(*)=OBSERVED, (X)=PREDICTED, (+)= IF * AND X OR'
!xx     A,' Y COINSIDE.')
!xx  933 FORMAT(1H ,50X,'(Y)=LONG RANGE FORECASTING')
!c  934 FORMAT(1H ,3X,E20.5,2(34X,E20.5))
!xx  934 FORMAT(/1H ,3X,D20.5,2(34X,D20.5))
!xx  935 FORMAT(1H ,10X,2H++,2(59(1H-),1H+))
!xx  936 FORMAT(1H ,' N =',I4,2X,1HI,121A1)
!c  937 FORMAT(1H ,14X,E20.5)
!xx  937 FORMAT(1H ,14X,D20.5)
!xx 9370 FORMAT(1H ,'   N=',I5,5(4X,D20.5))
!xx 9371 FORMAT(1H ,10X,5(4X,D20.5))
!xx 1357 FORMAT(/1H ,'MATRIX S')
!xx 8121 FORMAT(/1H ,'MATRIX  A(',I3,')')
return
end
!
!
subroutine maxmin(dmax,dmin,yd)
  use timsac_kinds, only: dp
  implicit none
!xx      REAL*8 DMAX,DMIN,YD
real(dp) dmax, dmin, yd
!c 1014 IF(DMAX-YD) 1001,1002,1002
!xx 1014 IF(DMAX-YD.GE.0) GO TO 1002
!xx 1001 DMAX=YD
if(dmax-yd.ge.0) go to 1002
dmax=yd
!c 1002 IF(DMIN-YD) 1004,1004,1003
1002 if(dmin-yd.le.0) go to 1004
!xx 1003 DMIN=YD
dmin=yd
1004 return
end
