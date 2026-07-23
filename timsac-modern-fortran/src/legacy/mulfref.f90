! Part of the experimental modern Fortran port of timsac 1.3.8-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original TIMSAC credits are retained; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

subroutine mulfrff(k,inw,n,lagh1,ip0,p,x,c,s,g,ph,pch,r,chm)
  use timsac_kinds, only: dp
  implicit none
!
!     PROGRAM 5.2.4   FREQUENCY RESPONSE FUNCTION (MULTIPLE CHANNEL)
!-----------------------------------------------------------------------
!      SUBROUTINE FQCPIV(X,XDET,MM,MJ)
!      SUBROUTINE MPHASE(C,S,OARC,PH,K,JJF)
!      SUBROUTINE MULARC(C,S,ARC,K)
!      SUBROUTINE MULERR(PCH,R,N,LAGH1,K,JJF,D1,D2)
!      SUBROUTINE MULPAC(ARC,OARC,PH,K,JJF)
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
!     THIS PROGRAM COMPUTES MULTIPLE FREQUENCY RESPONSE FUNCTION, GAIN,
!     PHASE, MULTIPLE COHERENCY, PARTIAL COHERENCY AND RELATIVE ERROR
!     STATISTICS.
!     A CARD WITH THE TOATL NUMBER(K) OF INPUT VARIABLES AND ANOTHER
!     WITH SPECIFICATION OF INPUT VARIABLES(INW(I),I=1,K) AND OUTPUT
!     VARIABLE(INW(K+1)) SHOULD BE ADDED ON TOP OF THE OUTPUT OF
!     PROGRAM 5.2.2 MULSPE TO FORM THE INPUT TO THIS PROGRAM.
!     WITHIN IP0 VARIABLES OF MULSPE OUTPUT, ONLY THOSE K+1 INW(I)-TH
!     VARIABLES ARE TAKEN INTO COMPUTATION.
!
!xx      IMPLICIT REAL*8(A-H,O-W)
!xx      IMPLICIT COMPLEX*16(X-Z)
!      DIMENSION P(10,10),X(10,10),C(10),S(10),G(10)
!      DIMENSION OARC(10),PH(10),PCH(10),R(10),INW(10)
!xx      DIMENSION P(LAGH1,IP0,IP0),X(IP0,IP0,LAGH1)
!xx      DIMENSION C(K,LAGH1),S(K,LAGH1),G(K,LAGH1)
!xx      DIMENSION OARC(K),PH(K,LAGH1),PCH(K,LAGH1),R(K,LAGH1),INW(K+1)
!xx      DIMENSION CHM(LAGH1)
!xx      DIMENSION XFR(IP0,IP0,LAGH1)
integer k, inw(k+1), n, lagh1, ip0
real(dp) p(lagh1,ip0,ip0), c(k,lagh1), s(k,lagh1),&
&g(k,lagh1), ph(k,lagh1), pch(k,lagh1),&
&r(k,lagh1), chm(lagh1)
complex(kind(0.0d0)) x(ip0,ip0,lagh1)
! local
integer i, ii, im1, j, jf, jjf, k1
real(dp) oarc(k), p00, ep, g2, g3, d1, d2
complex(kind(0.0d0)) xfr(ip0,ip0,lagh1), xdet
!
!     INPUT / OUTPUT DATA FILE OPEN
!      CALL SETWND
!      CALL FLOPN2(NFL)
!      IF (NFL.EQ.0) GO TO 999
!     ABSOLUTE DIMENSION USED FOR SUBROUTINE CALL
!      MJ=10
!     INPUT OUTPUT VARIABLE SPECIFICATION
!      READ(5,1) K
k1=k+1
!      READ(5,1) (INW(I),I=1,K1)
!     FOLLOWING INPUTS ARE OUTPUTS OF PROGRAM 5.2.2 MULSPE.
!      READ(5,1) N,LAGH,IP0
!      LAGH1=LAGH+1
!     INITIAL CONDITION PRINT OUT
!      WRITE(6,55)
!      WRITE(6,56)
!      WRITE(6,57) N,LAGH,K
!      WRITE(6,259) (INW(I),I=1,K1)
!     COMPUTATION START
do 10 jf=1,lagh1
jjf=jf
!      JFM1=JF-1
!      WRITE(6,58) JFM1
!     SPECTRUM INPUT
!      CALL REMATX(P,IP0,IP0,1,MJ,MJ)
!     REAL TO COMPLEX TRANSFORMATION
do 401 i=1,ip0
!      X(I,I)=P(I,I)
x(i,i,jf)=p(jf,i,i)
if(i.eq.1) go to 401
im1=i-1
do 402 j=1,im1
!      X(I,J)=DCMPLX(P(I,J),P(J,I))
!  402 X(J,I)=DCONJG(X(I,J))
!xx      X(I,J,JF)=DCMPLX(P(JF,I,J),P(JF,J,I))
x(i,j,jf)=cmplx(p(jf,i,j),p(jf,j,i),kind(0.0d0))
!xx  402 X(J,I,JF)=DCONJG(X(I,J,JF))
!xx      X(J,I,JF)=DCONJG(X(I,J,JF))
x(j,i,jf)=conjg(x(i,j,jf))
402 continue
401 continue
!     MATRIX REARRANGEMENT AND PRINT OUT (COMPLEX)
!      CALL REARRA(X,INW,IP0,K1,MJ)
call rearrac(x(1,1,jf),inw,ip0,k1)
!      WRITE(6,159)
!      CALL PRCPMA(X,K1,K1,MJ,MJ)
!     FREQUENCY RESPONSE FUNCTION COMPUTATION
!      P00=DREAL(X(K1,K1))
!      CALL FQCPIV(X,XDET,K,MJ)
!xx      P00=DREAL(X(K1,K1,JF))
p00=real(x(k1,k1,jf))
!      CALL FQCPIV(X(1,1,JF),XDET,K,MJ)
do 31 i=1,ip0
do 30 ii=1,ip0
xfr(i,ii,jf)=x(i,ii,jf)
30 continue
31 continue
call fqcpiv(xfr(1,1,jf),xdet,k,ip0)
do 20 i=1,k
!      C(I)=DREAL(X(I,K1))
!   20 S(I)=-DIMAG(X(I,K1))
!xx      C(I,JF)=DREAL(XFR(I,K1,JF))
c(i,jf)=real(xfr(i,k1,jf))
!xx   20 S(I,JF)=-DIMAG(XFR(I,K1,JF))
!xx      S(I,JF)=-DIMAG(XFR(I,K1,JF))
s(i,jf)=-aimag(xfr(i,k1,jf))
20 continue
!     GAIN COMPUTATION
do 21 i=1,k
!   21 G(I)=DSQRT(C(I)**2+S(I)**2)
!xx   21 G(I,JF)=DSQRT(C(I,JF)**2+S(I,JF)**2)
g(i,jf)=dsqrt(c(i,jf)**2+s(i,jf)**2)
21 continue
!     PHASE COMPUTATION
!      CALL MPHASE(C,S,OARC,PH,K,JJF)
if(jjf.ne.1) then
do 24 i=1,k
!xx   24    PH(I,JF)=PH(I,JF-1)
ph(i,jf)=ph(i,jf-1)
24 continue
end if
call mphase(c(1,jf),s(1,jf),oarc,ph(1,jf),k,jjf)
!     PARTIAL COHERENCY AND MULTIPLE COHERENCY COMPUTATION
!      EP=DREAL(X(K1,K1))
!xx      EP=DREAL(XFR(K1,K1,JF))
ep=real(xfr(k1,k1,jf))
do 22 i=1,k
!      G2=G(I)**2
!      G3=G2+EP*X(I,I)
g2=g(i,jf)**2
!xx      G3=G2+EP*XFR(I,I,JF)
!xx      G3=G2+EP*DREAL(XFR(I,I,JF))
g3=g2+ep*real(xfr(i,i,jf))
if(g3.ne.0.0) go to 23
!      PCH(I)=100.0D-00
pch(i,jf)=100.0d-00
go to 22
!   23 PCH(I)=G2/G3
23 pch(i,jf)=g2/g3
22 continue
!      CHM=1.0D-00-EP/P00
chm(jf)=1.0d-00-ep/p00
!     RELATIVE ERROR STATISTICS COMPUTATION
!      CALL MULERR(PCH,R,N,LAGH1,K,JJF,D1,D2)
call mulerr(pch(1,jf),r(1,jf),n,lagh1,k,jjf,d1,d2)
!     FREQUENCY RESPONSE FUNCTION, GAIN, PHASE, PARTIAL COHERENCY,
!     MULTIPLE COHERENCY, RELATIVE ERROR STATISTICS PRINT OUT
!      WRITE(6,60)
!      WRITE(6,61)
!      CALL PRCOL6(C,S,G,PH,PCH,R,1,K,0)
!      WRITE(6,65) CHM
!      WRITE(6,65) CHM(JF)
10 continue
!      CALL FLCLS2(NFL)
!  999 CONTINUE
!    1 FORMAT(10I5)
!   55 FORMAT(1H ,62HPROGRAM 5.2.4   FREQUENCY RESPONSE FUNCTION (MULTIPL
!     AE CHANNEL))
!   56 FORMAT(1H ,17HINITIAL CONDITION)
!   57 FORMAT(1H ,2HN=,I5,5X,5HLAGH=,I5,5X,2HK=,I5)
!   58 FORMAT(//1H ,2HF=,I5)
!   60 FORMAT(//1H ,4X,1HI,3X,27HFREQUENCY RESPONSE FUNCTION,10X,4HGAIN,9
!     AX,5HPHASE,7X,7HPARTIAL,6X,8HRELATIVE)
!   61 FORMAT(1H ,12X,9HREAL PART,4X,10HIMAG. PART,33X,9HCOHERENCY,9X,5HE
!     ARROR)
!   65 FORMAT(1H ,69X,8HMULTIPLE/1H ,68X,9HCOHERENCY/1H ,63X,D14.5)
!  159 FORMAT(1H ,28HSPECTRUM MATRIX (REARRANGED))
!  259 FORMAT(/1H ,6HINW(I),5X,10I5)
return
end subroutine
!
subroutine fqcpiv(x,xdet,mm,mj)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES MULTIPLE FREQUENCY RESPONSE FUNCTION.
!     MM: THE TOTAL NUMBER OF INPUTS (LESS THAN 10)
!     MJ: ABSOLUTE DIMENSION OF X IN THE MAIN ROUTINE
!xx      IMPLICIT COMPLEX*16(X)
!xx      DIMENSION X(MJ,MJ)
!x      DIMENSION IDS(10)
!xx      DIMENSION IDS(MM)
integer mm, mj
complex(kind(0.0d0)) x(mj,mj), xdet
! local
integer i, j, jj, l, maxi, mm1, mmj, mp1, ids(mm)
complex(kind(0.0d0)) xmaxp, xc
!
xdet=1.0d-00
mp1=mm+1
do 10 l=1,mm
!     PIVOTING AT L-TH STAGE
xmaxp=0.10000d-10
maxi=0
do 110 i=l,mm
!xx      IF(CDABS(XMAXP).GE.CDABS(X(I,L))) GO TO 110
if(abs(xmaxp).ge.abs(x(i,l))) go to 110
xmaxp=x(i,l)
maxi=i
110 continue
ids(l)=maxi
if(maxi.eq.l) go to 120
if(maxi.gt.0) go to 121
xdet=0.0d-00
go to 140
!     ROW INTERCHANGE
121 do 14 j=1,mp1
xc=x(maxi,j)
x(maxi,j)=x(l,j)
!xx   14 X(L,J)=XC
x(l,j)=xc
14 continue
xdet=-xdet
120 xdet=xdet*xmaxp
xc=1.0d-00/xmaxp
x(l,l)=1.0d-00
do 11 j=1,mp1
!xx   11 X(L,J)=X(L,J)*XC
x(l,j)=x(l,j)*xc
11 continue
do 12 i=1,mp1
if(i.eq.l) go to 12
xc=x(i,l)
x(i,l)=0.0d-00
do 13 j=1,mp1
!xx   13 X(I,J)=X(I,J)-XC*X(L,J)
x(i,j)=x(i,j)-xc*x(l,j)
13 continue
12 continue
10 continue
if(mm.gt.1) go to 123
go to 140
!     COLUMN INTERCHANGE
123 mm1=mm-1
do 130 j=1,mm1
mmj=mm-j
jj=ids(mmj)
if(jj.eq.mmj) go to 130
do 131 i=1,mp1
xc=x(i,jj)
x(i,jj)=x(i,mmj)
!xx  131 X(I,MMJ)=XC
x(i,mmj)=xc
131 continue
130 continue
140 return
end subroutine
!
subroutine mphase(c,s,oarc,ph,k,jjf)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES PHASE.
!     (MULTIPLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION C(K),S(K),OARC(K),PH(K)
!      DIMENSION ARC(10)
!xx      DIMENSION ARC(K)
integer k, jjf
real(dp) c(k), s(k), oarc(k), ph(k)
real(dp) arc(k)
!     ARCTANGENT COMPUTATION
call mularc(c,s,arc,k)
!     PHASE COMPUTATION
call mulpac(arc,oarc,ph,k,jjf)
return
end subroutine
!
subroutine mularc(c,s,arc,k)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES RAW PHASE.
!     (MULTIPLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION C(K),S(K),ARC(K)
integer k
real(dp) c(k), s(k), arc(k)
! local
integer i
real(dp) pi, cst5
pi=3.1415926536
cst5=0.5d-00
do 10 i=1,k
!c      IF(C(I)) 11,12,13
!c   11 IF(S(I)) 14,15,16
!c   12 IF(S(I)) 17,18,19
if(c(i).eq.0) go to 12
if(c(i).gt.0) go to 13
!xx   11 IF(S(I).LT.0) GO TO 14
if(s(i).lt.0) go to 14
if(s(i).eq.0) go to 15
if(s(i).gt.0) go to 16
12 if(s(i).lt.0) go to 17
if(s(i).eq.0) go to 18
if(s(i).gt.0) go to 19
13 arc(i)=datan(s(i)/c(i))
go to 10
14 arc(i)=datan(s(i)/c(i))-pi
go to 10
15 arc(i)=-pi
go to 10
16 arc(i)=datan(s(i)/c(i))+pi
go to 10
17 arc(i)=-pi*cst5
go to 10
18 arc(i)=0.0d-00
go to 10
19 arc(i)=pi*cst5
10 continue
return
end subroutine
!
subroutine mulerr(pch,r,n,lagh1,k,jjf,d1,d2)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE COMPUTES RELATIVE ERROR STATISTICS.
!     (MULTIPLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION PCH(K),R(K)
integer n, lagh1, k,jjf
real(dp) pch(k), r(k), d1, d2
! local
integer i, lagh
real(dp) cst0, cst1, cst100, e1, er
cst0=0.0d-00
cst1=1.0d-00
cst100=100.0d-00
if(jjf.ne.1) go to 30
!     CONSTANTS D1,D2 COMPUTATION
lagh=lagh1-1
call subd12(n,lagh,k,d1,d2)
!     RELATIVE ERROR STATISTICS COMPUTATION
30 do 20 i=1,k
if(pch(i).le.cst0) go to 22
if(pch(i).gt.cst1) go to 22
e1=cst1/pch(i)-cst1
er=dsqrt(e1)
if(jjf.eq.1) go to 23
if(jjf.eq.lagh1) go to 23
r(i)=d2*er
go to 20
23 r(i)=d1*er
go to 20
22 r(i)=cst100
20 continue
return
end subroutine
!
subroutine mulpac(arc,oarc,ph,k,jjf)
  use timsac_kinds, only: dp
  implicit none
!     THIS SUBROUTINE MAKES PHASE CURVE CONTINUOUS.
!     (MULTIPLE CHANNEL)
!xx      IMPLICIT REAL*8(A-H,O-Z)
!xx      DIMENSION ARC(K),OARC(K),PH(K)
integer k, jjf
real(dp) arc(k), oarc(k), ph(k)
! local
integer i
real(dp) pi, pi2, dk
pi=3.1415926536
pi2=pi+pi
if(jjf.ne.1) go to 20
do 9 i=1,k
ph(i)=arc(i)
!xx    9 OARC(I)=ARC(I)
oarc(i)=arc(i)
9 continue
go to 30
!xx   20 DO 10 I=1,K
20 do 100 i=1,k
dk=arc(i)-oarc(i)
if(dk.gt.pi) go to 11
if(dk.lt.-pi) go to 12
ph(i)=ph(i)+dk
go to 10
11 ph(i)=ph(i)+dk-pi2
go to 10
12 ph(i)=ph(i)+dk+pi2
10 oarc(i)=arc(i)
100 continue
30 continue
return
end subroutine
