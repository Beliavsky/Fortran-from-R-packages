! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mulrspf(h,l,ip,k,sd,a,b,y,ch)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.4.2   MULTIPLE RATIONAL SPECTRUM
!-----------------------------------------------------------------------
!      SUBROUTINE XYCTRX(X,Y,Z,MM,NN)
!-----------------------------------------------------------------------
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: MARCH 25, 1977
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
!         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
!         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
!-----------------------------------------------------------------------
!     THIS PROGRAM COMPUTES RATIONAL SPECTRUM FOR IP-DIMENSIONAL
!     AR-MA PROCESS
!     X(N)=A(1)X(N-1)+...+A(L)X(N-L)+E(N)+B(1)E(N-1)+...+B(K)E(N-K),
!     WHERE E(N) IS A WHITE NOISE WITH ZERO MEAN VECTOR AND COVARIANCE
!     MATRIX SD.
!     OUTPUTS ARE SPECTRUM MATRIX P(I) AT FREQUENCIES I/(2*H)
!     (I=0,1,...,H).
!
!xx      IMPLICIT REAL*8(A-H,O-W)
!xx      IMPLICIT COMPLEX*16(X-Z)
!xx      INTEGER H,H1
!      COMMON G,GR,GI,LG,H,JJF
!      DIMENSION SD(10,10),A(30,10,10),B(30,10,10),X(10,10),Y(10,10)
!      DIMENSION Z(10,10),G(31),CH(10,10)
!xx      DIMENSION SD(IP,IP),A(L,IP,IP),B(K,IP,IP),X(IP,IP),Y(IP,IP,H+1)
!xx      DIMENSION Z(IP,IP),G(L+K+1),CH(IP,IP,H+1)
integer h, l, ip, k
real(dp) sd(ip,ip), a(l,ip,ip), b(k,ip,ip), ch(ip,ip,h+1)
complex(kind(0.0d0)) y(ip,ip,h+1)
! local
integer i, i1, ii, im1, jf, jj, jjf, lg, h1
real(dp) g(l+k+1), cst0,  cst1, gr, gi, ryi, ryj, rryij,&
&riyij
complex(kind(0.0d0)) x(ip,ip), z(ip,ip), xdet
!     INPUT / OUTPUT DATA FILE OPEN
!	CALL SETWND
!	CALL FLOPN2(NFL)
!	IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSIONS USED FOR SUBROUTINE CALL
!      MJ0=30
!      MJ1=10
cst0=0.0d-00
cst1=1.0d-00
!     H SPECIFICATION
!      READ(5,1) H
!     SD AND A INPUT
!     THE OUTPUTS OF PROGRAM 5.3.2 FPEC(WITH IL=0) CAN BE USED AS THE
!     FOLLWOING INPUTS WITH K=0.
!      READ(5,1) N,L,IP
!     SD INPUT
!      CALL REMATX(SD,IP,IP,1,MJ1,MJ1)
!      IF(L.LE.0) GO TO 300
!     A INPUT
!      CALL REMAT3(A,L,IP,IP,1,MJ0,MJ1,MJ1)
!     K INPUT
!  300 READ(5,1) K
!      IF(K.LE.0) GO TO 310
!     B INPUT
!      CALL REMAT3(B,K,IP,IP,1,MJ0,MJ1,MJ1)
!xx  310 H1=H+1
h1=h+1
!     INITIAL CONDITION PRINT OUT
!      WRITE(6,59)
!      WRITE(6,60)
!      WRITE(6,61) H,N,L,IP,K
!      WRITE(6,161)
!      CALL SUBMPR(SD,IP,IP,MJ1,MJ1)
!      IF(L.LE.0) GO TO 400
!     A PRINT OUT
!      WRITE(6,420)
!      CALL PRMAT3(A,L,IP,IP,0,MJ0,MJ1,MJ1)
!  400 IF(K.LE.0) GO TO 410
!     B PRINT OUT
!      WRITE(6,430)
!      CALL PRMAT3(B,K,IP,IP,0,MJ0,MJ1,MJ1)
!     SPECTRUM COMPUTATION
!xx  410 DO 10 JF=1,H1
do 10 jf=1,h1
jjf=jf
!     SD STORE
!xx      DO 631 II=1,IP
do 632 ii=1,ip
do 631 jj=1,ip
!  631 Y(II,JJ)=SD(II,JJ)
!xx  631 Y(II,JJ,JF)=SD(II,JJ)
y(ii,jj,jf)=sd(ii,jj)
631 continue
632 continue
if(k.gt.0) go to 100
!xx      DO 110 II=1,IP
do 111 ii=1,ip
do 110 jj=1,ip
!xx  110 Z(II,JJ)=SD(II,JJ)
z(ii,jj)=sd(ii,jj)
110 continue
111 continue
go to 224
!     BF COMPUTATION
100 do 20 ii=1,ip
do 21 jj=1,ip
if(ii.ne.jj) go to 22
g(1)=cst1
go to 23
22 g(1)=cst0
23 do 25 i=1,k
i1=i+1
!xx   25 G(I1)=B(I,II,JJ)
g(i1)=b(i,ii,jj)
25 continue
!xx   24 LG=K
lg=k
!      CALL FGER1
call fger1(g,gr,gi,lg,h,jjf)
!xx      X(II,JJ)=DCMPLX(GR,GI)
x(ii,jj)=cmplx(gr,gi,kind(0.0d0))
21 continue
20 continue
!     BF*SD*CONJG(BF') COMPUTATION
!      CALL XYCTRX(X,Y,Z,IP,IP,MJ1,MJ1)
call xyctrx(x,y(1,1,jf),z,ip,ip)
224 if(l.gt.0) go to 120
!xx      DO 130 II=1,IP
do 131 ii=1,ip
do 130 jj=1,ip
!  130 Y(II,JJ)=Z(II,JJ)
!xx  130 Y(II,JJ,JF)=Z(II,JJ)
y(ii,jj,jf)=z(ii,jj)
130 continue
131 continue
go to 244
!     AF COMPUTATION
120 do 40 ii=1,ip
do 41 jj=1,ip
if(ii.ne.jj) go to 42
g(1)=cst1
go to 43
42 g(1)=cst0
43 do 45 i=1,l
i1=i+1
!xx   45 G(I1)=-A(I,II,JJ)
g(i1)=-a(i,ii,jj)
45 continue
!xx   44 LG=L
lg=l
!      CALL FGER1
call fger1(g,gr,gi,lg,h,jjf)
!xx      X(II,JJ)=DCMPLX(GR,GI)
x(ii,jj)=cmplx(gr,gi,kind(0.0d0))
41 continue
40 continue
!     INVERSE OF AF (COMPLEX) COMPUTATION
!      CALL INVDET(X,XDET,IP,MJ1)
call invdetc(x,xdet,ip)
!     (INVERSE OF AF)*(BF*SD*CONJG(BF'))*CONJG((INVERSE OF AF)')
!     COMPUTATION
!      CALL XYCTRX(X,Z,Y,IP,IP,MJ1,MJ1)
call xyctrx(x,z,y(1,1,jf),ip,ip)
!     SIMPLE COHERENCE COMPUTATION
!  244 CH(1,1)=CST1
244 ch(1,1,jf)=cst1
!c      IF(IP.EQ.1) GO TO 260
if(ip.eq.1) go to 10
do 50 ii=2,ip
im1=ii-1
!      RYI=DREAL(Y(II,II))
!xx      RYI=DREAL(Y(II,II,JF))
ryi=real(y(ii,ii,jf))
do 51 jj=1,im1
!      RYJ=DREAL(Y(JJ,JJ))
!      RRYIJ=DREAL(Y(II,JJ))
!      RIYIJ=DIMAG(Y(II,JJ))
!      CH(II,JJ)=(RRYIJ**2+RIYIJ**2)/(RYI*RYJ)
!   51 CH(JJ,II)=CH(II,JJ)
!   50 CH(II,II)=CST1
!xx      RYJ=DREAL(Y(JJ,JJ,JF))
!xx      RRYIJ=DREAL(Y(II,JJ,JF))
!xx      RIYIJ=DIMAG(Y(II,JJ,JF))
ryj=real(y(jj,jj,jf))
rryij=real(y(ii,jj,jf))
riyij=aimag(y(ii,jj,jf))
ch(ii,jj,jf)=(rryij**2+riyij**2)/(ryi*ryj)
!xx   51 CH(JJ,II,JF)=CH(II,JJ,JF)
!xx   50 CH(II,II,JF)=CST1
ch(jj,ii,jf)=ch(ii,jj,jf)
51 continue
ch(ii,ii,jf)=cst1
50 continue
!     RATIONAL SPECTRUM AND SIMPLE COHERENCE PRINT OUT
!c  260 JFM1=JF-1
!      WRITE(6,65) JFM1
!      WRITE(6,66)
!      CALL PRCPMA(Y,IP,IP,MJ1,MJ1)
!      WRITE(6,67)
!      CALL SUBMPR(CH,IP,IP,MJ1,MJ1)
!
10 continue
!	CALL FLCLS2(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!   59 FORMAT(1H ,42HPROGRAM 5.4.2   MULTIPLE RATIONAL SPECTRUM)
!   60 FORMAT(1H ,17HINITIAL CONDITION)
!   61 FORMAT(1H ,2HH=,I5,5X,2HN=,I5,5X,2HL=,I5,5X,3HIP=,I5,5X,2HK=,I5)
!   62 FORMAT(1H ,2HI=,I5)
!xx   65 FORMAT(///1H ,2HF=,I5)
!xx   66 FORMAT(1H ,5X,17HRATIONAL SPECTRUM)
!xx   67 FORMAT(/1H ,5X,16HSIMPLE COHERENCE)
!  161 FORMAT(//1H ,7HSD(I,J))
!  420 FORMAT(//1H ,6HA(I,J))
!  430 FORMAT(//1H ,6HB(I,J))
!    2 FORMAT(4D20.10)
return
end
!
!      SUBROUTINE XYCTRX(X,Y,Z,MM,NN,MJ1,MJ2)
subroutine xyctrx(x,y,z,mm,nn)
  use timsac_kinds, only: dp
  implicit none
!     Z=X*Y*CONJG(X')
!     Y,Z: HERMITIAN
!     (UPPER LEFT MM X MM OF Z)=(UPPER LEFT MM X NN OF X)*(UPPER LEFT
!     NN X NN OF Y)*CONJG((UPPER LEFT MM X NN OF X)')
!     (MJ1,MJ2): ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!     (MJ2,MJ2): ABSOLUTE DIMENSION OF Y IN THE MAIN ROUTINE
!     (MJ1,MJ1): ABSOLUTE DIMENSION OF Z IN THE MAIN ROUTINE
!     MM,NN: SHOULD BE LESS THAN 11.
!xx      IMPLICIT COMPLEX*16(X-Z)
!      DIMENSION X(MJ1,MJ2),Y(MJ2,MJ2),Z(MJ1,MJ1)
!      DIMENSION Y1(10,10)
!xx      DIMENSION X(MM,NN),Y(NN,NN),Z(MM,MM)
!xx      DIMENSION Y1(MM,NN)
!xx      DOUBLE PRECISION CST0
integer mm, nn
complex(kind(0.0d0)) x(mm,nn), y(nn,nn), z(mm,mm)
! local
integer i, j, k
real(dp) cst0
complex(kind(0.0d0)) y1(mm,nn), xsum
cst0=0.0d-00
do 10 i=1,mm
!xx      DO 10 J=1,NN
do 11 j=1,nn
xsum=cst0
do 12 k=1,nn
!xx   12 XSUM=XSUM+X(I,K)*Y(K,J)
xsum=xsum+x(i,k)*y(k,j)
12 continue
!xx   10 Y1(I,J)=XSUM
y1(i,j)=xsum
11 continue
10 continue
do 110 i=1,mm
!xx      DO 110 J=1,I
do 111 j=1,i
xsum=cst0
do 112 k=1,nn
!xx  112 XSUM=XSUM+Y1(I,K)*DCONJG(X(J,K))
!xx      XSUM=XSUM+Y1(I,K)*DCONJG(X(J,K))
xsum=xsum+y1(i,k)*conjg(x(j,k))
112 continue
z(i,j)=xsum
!xx  110 Z(J,I)=DCONJG(Z(I,J))
!xx      Z(J,I)=DCONJG(Z(I,J))
z(j,i)=conjg(z(i,j))
111 continue
110 continue
return
end
