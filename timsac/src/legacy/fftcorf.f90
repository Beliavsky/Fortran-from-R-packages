! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine fftcorf(ld,lagh1,n,n2p,isw,x1,y1,xa,x,y,&
&cna1,cn1,cn2,amean)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.1.3   AUTO AND/OR CROSS CORRELATIONS VIA FFT.
!-----------------------------------------------------------------------
!     ** THIS PROGRAM IS AN ADAPTED VERSION OF THE ALGOL PROCEDURE
!         //CROSS CORRELATE// PREPARED BY DR GORDON SANDE AT PRINCETON
!         UNIVERSITY IN 1966-7.
!     ** DESIGNED BY H. AKAIKE, THE INSTITUTE OF STATISTICAL MATHEMATICS
!     ** PROGRAMMED BY E. ARAHATA, THE INSTITUTE OF STATISTICAL MATHEMAT
!         TOKYO
!     ** DATE OF THE LATEST REVISION: FEB. 16, 1978
!     ** THIS PROGRAM WAS ORIGINALLY PUBLISHED IN
!         "DAINAMIKKU SISTEMU NO TOKEI-TEKI KAISEKI TO SEIGYO (STATISTICA
!         ANALYSIS AND CONTROL OF DYNAMIC SYSTEMS)" BY H. AKAIKE AND
!         T. NAKAGAWA, SAIENSU-SHA, TOKYO, 1972 (IN JAPANESE)
!-----------------------------------------------------------------------
!     THIS PROGRAM COMPUTES AUTO AND/OR CROSS
!     COVARIANCES AND CORRELATIONS VIA FFT.
!     IT REQUIRES FOLLOWING INPUTS:
!     ISW: ISW=1...AUTO CORRELATION OF X (ONE-CHANNEL)
!           ISW=2...AUTO CORRELATIONS OF X AND Y (TWO-CHANNEL)
!           ISW=4...AUTO,CROSS CORRELATIONS OF X AND Y (TWO-CHANNEL)
!     LD: LENGTH OF DATA
!     LAGH: MAXIMUM LAG
!     DFORM: INPUT FORMAT SPECIFICATION STATEMENT IN ONE CARD,
!     FOR EXAMPLE
!     (8F10.4)
!     (X(I); I=1,LD): DATA OF CHANNEL X
!     (Y(I); I=1,LD): DATA OF CHANNEL Y (FOR ISW=2 OR 4 ONLY)
!
!xx      IMPLICIT REAL*8(A-H,O-Y)
!xx      IMPLICIT COMPLEX*16(Z)
!      REAL*4 DFORM
!      DIMENSION X(2048),Y(2048),Z(2048),ZS(1025)
!      DIMENSION CN1(501),CN2(501)
!      DIMENSION DFORM(20)
!      REAL*4 XS,YS
!      DIMENSION XS(2048),YS(2048)
!xx      DIMENSION X1(LD),Y1(LD),XA(N,2),X(N),Y(N)
!xx      DIMENSION Z(N),ZS(N/2+1)
!xx      DIMENSION CNA1(LAGH1,2),CN1(LAGH1),CN2(LAGH1),AMEAN(2)
integer ld, lagh1, n, n2p, isw
real(dp) x1(ld), y1(ld), xa(n,2), x(n), y(n),&
&cna1(lagh1,2), cn1(lagh1), cn2(lagh1), amean(2)
! local
integer i, i1, ii, isg, j1, lagh, m, m1, nd, ni, np1, np2
real(dp) cst0, cst1, cst2, an, ald, ald1, rf, sf, rg,&
&sg, xi, xni, yi, yni, x0, xmean, ymean, cx0, y0
complex(kind(0.0d0)) z(n), zs(n/2+1), zi, zni
!     INPUT / OUTPUT DATA FILE OPEN
!	CHARACTER(100) DFNAM
!	DFNAM='fftcor.out'
!	CALL SETWND
!	CALL FLOPN3(DFNAM,NFL)
!	IF (NFL.EQ.0) GO TO 999
!     INITIAL CONDITION INPUT AND PUNCH OUT
!      READ(5,1) ISW,LD,LAGH
!      WRITE(6,50)
!      WRITE(6,51)
!      WRITE(6,52) ISW,LD,LAGH
!      WRITE(7,1) LD,LAGH
!      LAGH1=LAGH+1
lagh=lagh1-1
nd=ld+lagh1
!     N2P, N: DEFINITION
!      I0=1
!   10 IR1=2**I0
!      IF(IR1-ND) 11,12,12
!   11 I0=I0+1
!      GO TO 10
!   12 N2P=I0
!      N=2**N2P
np1=n+1
np2=n+2
m=n/2
m1=m+1
cst0=0.0d-00
cst1=1.0d-00
cst2=0.25d-00
an=n
ald=ld
ald1=cst1/(an*ald)
!     INPUT FORMAT SPECIFICATION
!      READ(5,4) (DFORM(I),I=1,20)
!    4 FORMAT(20A4)
!     ORIGINAL DATA INPUT AND OUTPUT
do 20 i=1,n
x(i)=cst0
y(i)=cst0
20 continue
!      READ(5,DFORM) (XS(I),I=1,LD)
do 1200 i=1,ld
!      X(I)=DBLE(XS(I))
x(i)=x1(i)
1200 continue
if(isw.eq.1) go to 200
!      READ(5,DFORM) (YS(I),I=1,LD)
do 1201 i=1,ld
!      Y(I)=DBLE(YS(I))
y(i)=y1(i)
1201 continue
!  200 WRITE(6,53)
200 continue
!      IF(ISW.NE.1) GO TO 201
!      WRITE(6,54)
!      CALL PRCOL1(X,1,LD,0)
!      GO TO 202
!  201 WRITE(6,55)
!      CALL PRCOL2(X,Y,1,LD,0)
!     MEAN DELETION
!xx  202 CALL DMEADL(X,LD,XMEAN)
call dmeadl(x,ld,xmean)
if(isw.eq.1) go to 203
call dmeadl(y,ld,ymean)
!     DOUBLE PRECISION COMPLEX REPRESENTATION
203 do 31 i=1,n
!xx   31 Z(I)=DCMPLX(X(I),Y(I))
!xx      Z(I)=DCMPLX(X(I),Y(I))
z(i)=cmplx(x(i),y(i),kind(0.0d0))
31 continue
!     FOURIER TRANSFORM OF Z
isg=-1
call mixrad(z,n,n2p,isg)
if(isw.ne.1) go to 204
!     RAW SPECTRUM COMPUTATION
do 32 i=2,m
!xx      X(I)=DREAL(Z(I))**2+DIMAG(Z(I))**2
x(i)=real(z(i))**2+aimag(z(i))**2
ni=np2-i
!xx   32 X(NI)=X(I)
x(ni)=x(i)
32 continue
!xx      X(1)=DREAL(Z(1))**2
!xx      X(M1)=DREAL(Z(M1))**2
x(1)=real(z(1))**2
x(m1)=real(z(m1))**2
go to 205
!     DECOMPOSITION AND RAW SPECTRUM COMPUTATION
204 do 125 i=2,m
ni=np2-i
zi=z(i)
zni=z(ni)
!xx      RF=DREAL(ZI)
!xx      SF=DIMAG(ZI)
!xx      RG=DREAL(ZNI)
!xx      SG=DIMAG(ZNI)
rf=real(zi)
sf=aimag(zi)
rg=real(zni)
sg=aimag(zni)
xi=rf+rg
xni=sf-sg
!xx      Z(I)=DCMPLX(XI,XNI)
z(i)=cmplx(xi,xni,kind(0.0d0))
x(i)=cst2*(xi**2+xni**2)
x(ni)=x(i)
yi=sf+sg
yni=rf-rg
!xx      Z(NI)=DCMPLX(YI,YNI)
z(ni)=cmplx(yi,yni,kind(0.0d0))
y(i)=cst2*(yi**2+yni**2)
y(ni)=y(i)
125 continue
!xx      X(1)=DREAL(Z(1))**2
!xx      Y(1)=DIMAG(Z(1))**2
!xx      X(M1)=DREAL(Z(M1))**2
!xx      Y(M1)=DIMAG(Z(M1))**2
x(1)=real(z(1))**2
y(1)=aimag(z(1))**2
x(m1)=real(z(m1))**2
y(m1)=aimag(z(m1))**2
if(isw.ne.4) go to 205
!     RAW CROSS SPECTRUM COMPUTATION
do 126 i=2,m
ni=np2-i
!xx  126 ZS(I)=CST2*Z(I)*Z(NI)
zs(i)=cst2*z(i)*z(ni)
126 continue
!xx      ZS(1)=DREAL(Z(1))*DIMAG(Z(1))
!xx      ZS(M1)=DREAL(Z(M1))*DIMAG(Z(M1))
zs(1)=real(z(1))*aimag(z(1))
zs(m1)=real(z(m1))*aimag(z(m1))
!     AUTO COVARIANCE COMPUTATION
205 do 33 i=1,n
!xx   33 Z(I)=DCMPLX(X(I),Y(I))
!xx      Z(I)=DCMPLX(X(I),Y(I))
z(i)=cmplx(x(i),y(i),kind(0.0d0))
33 continue
!     FOURIER TRANSFORM
!xx  215 CALL MIXRAD(Z,N,N2P,ISG)
call mixrad(z,n,n2p,isg)
ii=1
do 34 i=1,lagh1
!xx      X(I)=DREAL(Z(I))*ALD1
x(i)=real(z(i))*ald1
!xx   34 XA(I,II)=X(I)
xa(i,ii)=x(i)
34 continue
x0=x(1)
!      AMEAN=XMEAN
amean(ii)=xmean
!     NORMALIZATION
36 cx0=x(1)
!      CALL CORNOM(X,CN1,LAGH1,CX0,CX0)
call cornom(x,cna1(1,ii),lagh1,cx0,cx0)
!     AUTO COVARIANCE PRINT OUT
!      WRITE(6,162) II,II,AMEAN
!      WRITE(6,163)
!      CALL PRCOL2(X,CN1,1,LAGH1,1)
!     AUTO COVARIANCE PUNCH OUT
!      WRITE(7,1) II,II
!      WRITE(7,2) (X(I),I=1,LAGH1)
if(isw.eq.1) go to 300
if(ii.eq.2) go to 216
ii=2
do 35 i=1,lagh1
!xx      X(I)=DIMAG(Z(I))*ALD1
x(i)=aimag(z(i))*ald1
!xx   35 XA(I,II)=X(I)
xa(i,ii)=x(i)
35 continue
y0=x(1)
!      AMEAN=YMEAN
amean(ii)=ymean
go to 36
216 if(isw.ne.4) go to 300
!     CROSS COVARIANCE COMPUTATION
do 127 i=2,m
ni=np2-i
z(i)=zs(i)
!xx  127 Z(NI)=DCONJG(ZS(I))
!xx      Z(NI)=DCONJG(ZS(I))
z(ni)=conjg(zs(i))
127 continue
z(1)=zs(1)
z(m1)=zs(m1)
!     FOURIER TRANSFORM
call mixrad(z,n,n2p,isg)
do 41 i=1,lagh
i1=i+1
j1=np1-i
!xx      X(I1)=DREAL(Z(I1))*ALD1
x(i1)=real(z(i1))*ald1
!xx   41 Y(I1)=DREAL(Z(J1))*ALD1
!xx      Y(I1)=DREAL(Z(J1))*ALD1
y(i1)=real(z(j1))*ald1
41 continue
!xx      X(1)=DREAL(Z(1))*ALD1
x(1)=real(z(1))*ald1
y(1)=x(1)
!     NORMALIZATION
call cornom(x,cn1,lagh1,x0,y0)
call cornom(y,cn2,lagh1,x0,y0)
!     CROSS COVARIANCE PRINT OUT
!      JJ=1
!      WRITE(6,165) II,JJ
!      WRITE(6,166)
!      CALL PRCOL4(X,CN1,Y,CN2,1,LAGH1,1)
!     CROSS COVARIANCE PUNCH OUT
!      WRITE(7,1) II,JJ
!      WRITE(7,2) (X(I),I=1,LAGH1)
!      WRITE(7,1) JJ,II
!      WRITE(7,2) (Y(I),I=1,LAGH1)
300 continue
!	CALL FLCLS3(NFL)
!  999 CONTINUE
return
!    1 FORMAT(10I5)
!    2 FORMAT(4D20.10)
!   50 FORMAT(1H ,71HPROGRAM 5.1.3   AUTO AND/OR CROSS COVARIANCES AND CO
!     ARRELATIONS VIA FFT.)
!   51 FORMAT(1H ,17HINITIAL CONDITION)
!   52 FORMAT(1H ,4HISW=,I5,5X,3HLD=,I5,5X,5HLAGH=,I5)
!   53 FORMAT(1H ,13HORIGIANL DATA)
!   54 FORMAT(1H ,4X,1HI,12X,4HX(I))
!   55 FORMAT(1H ,4X,1HI,12X,4HX(I),10X,4HY(I))
!  162 FORMAT(//1H ,14HAUTOCOVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5,5
!     AX,5HMEAN=,D15.5)
!  163 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED)
!  165 FORMAT(//1H ,16HCROSS COVARIANCE,5X,6HCIJ(L),5X,2HI=,I5,5X,2HJ=,I5
!     A)
!  166 FORMAT(1H ,4X,1HL,5X,6HCIJ(L),8X,10HNORMALIZED,4X,6HCJI(L),8X,10HN
!     AORMALIZED)
end subroutine
